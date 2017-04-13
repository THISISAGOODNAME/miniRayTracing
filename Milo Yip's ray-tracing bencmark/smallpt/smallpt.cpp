#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <time.h>		

#ifdef __MSC_VER
#define INLINE __forceinline 
//#define INLINE __declspec(noinline)
#else
#define INLINE inline
#endif

#ifndef M_PI
#define M_PI 3.141592653589793238462643	
#endif

class RandomLCG {
	unsigned mSeed;

public:
	INLINE RandomLCG(unsigned seed = 0) : mSeed(seed) {
	}

	INLINE double operator()() { 
		mSeed = 214013u * mSeed + 2531011u; 
		return mSeed * (1.0 / 4294967296.0);
	}
};

struct Vec { 
	double x, y, z;

	static const Vec Zero;
	static const Vec XAxis;
	static const Vec YAxis;
	static const Vec ZAxis;

	INLINE Vec(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {
	}

	INLINE Vec operator+(const Vec& b) const { 
		return Vec(x + b.x, y + b.y,z + b.z); 
	}

	INLINE Vec operator-(const Vec& b) const { 
		return Vec(x - b.x, y - b.y, z - b.z); 
	}

	INLINE Vec operator*(double b) const { 
		return Vec(x * b, y * b, z * b); 
	}

	friend INLINE Vec operator*(const Vec& a, const Vec& b) {
		return Vec(a.x * b.x , a.y * b. y, a.z * b.z); 
	}

	INLINE Vec norm() const { 
		return *this * (1.0 / sqrt(x * x + y * y + z * z)); 
	}

	INLINE double operator%(const Vec& b) const { // dot
		return x * b.x + y * b.y + z * b.z; 
	}

	INLINE Vec operator^(const Vec& b) const { // cross
		return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
	}
};

const Vec Vec::Zero(0, 0, 0);
const Vec Vec::XAxis(1, 0, 0);
const Vec Vec::YAxis(0, 1, 0);
const Vec Vec::ZAxis(0, 0, 1);

struct Ray { 
	Vec o, d;

	INLINE Ray(const Vec &o, const Vec &d) : o(o), d(d) {
	} 
};

// material types, used in radiance()
enum Refl_t { 
	DIFF, 
	SPEC, 
	REFR 
};

struct Sphere {
	Vec p, e, c;      // position, emission, color
	Vec cc;
	double rad;       // radius
	double sqRad; 
	double maxC;
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)

	Sphere(double rad, const Vec& p, const Vec& e, const Vec& c, Refl_t refl) : p(p), e(e), c(c), rad(rad), refl(refl) {
		sqRad = rad * rad;
		maxC = c.x > c.y && c.y > c.z ? c.x : c.y > c.z ? c.y : c.z;
		cc = c * (1.0 / maxC);
	}

	// returns distance, 1e20 if nohit
	INLINE double intersect(const Ray &r) const	{
		// Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		Vec op = p - r.o;
		double b = op % r.d;
		double det = b * b - op % op + sqRad;
		const double eps = 1e-4;

		if (det < 0)
			return 1e20;
		else {
			double dets = sqrt(det);

			if (b - dets > eps)
				return b - dets;
			else if (b + dets > eps)
				return b + dets;
			else
				return 1e20;
		}
	}
};

//Scene: radius, position, emission, color, material
static Sphere spheres[] = {
	Sphere(1e5,  Vec( 1e5+1,40.8,81.6),  Vec::Zero, Vec(.75,.25,.25), DIFF),//Left
	Sphere(1e5,  Vec(-1e5+99,40.8,81.6), Vec::Zero, Vec(.25,.25,.75), DIFF),//Rght
	Sphere(1e5,  Vec(50,40.8, 1e5),      Vec::Zero, Vec(.75,.75,.75), DIFF),//Back
	Sphere(1e5,  Vec(50,40.8,-1e5+170),  Vec::Zero, Vec::Zero,        DIFF),//Frnt
	Sphere(1e5,  Vec(50, 1e5, 81.6),     Vec::Zero, Vec(.75,.75,.75), DIFF),//Botm
	Sphere(1e5,  Vec(50,-1e5+81.6,81.6), Vec::Zero, Vec(.75,.75,.75), DIFF),//Top
	Sphere(16.5, Vec(27,16.5,47),        Vec::Zero, Vec(1,1,1)*.999,  SPEC),//Mirr
	Sphere(16.5, Vec(73,16.5,78),        Vec::Zero, Vec(1,1,1)*.999,  REFR),//Glas
	Sphere(600,  Vec(50,681.6-.27,81.6), Vec(12,12,12), Vec::Zero,    DIFF) //Lite
};

INLINE double clamp(double x) { 
	if (x < 0)
		return 0;
	else if (x > 1)
		return 1;
	else
		return x;
}

INLINE int toInt(double x) { 
	return int(pow(clamp(x), 1 / 2.2) * 255 + .5); 
}

INLINE Sphere* intersect(const Ray &r, double &t) {
	t = 1e20;
	Sphere* ret = NULL;

	for (Sphere* s = spheres; s != spheres + sizeof(spheres) / sizeof(Sphere); ++s) {
		double d = s->intersect(r);
		if (d < t) {
			t = d;
			ret = s;
		}
	}
	return ret;
}

static Vec radiance(const Ray &r, int depth, RandomLCG& rand) {
	double t;                               // distance to intersection
	Sphere* obj = intersect(r, t);

	if (!obj)
		return Vec::Zero; // if miss, return black
	else {
		int newDepth = depth + 1;
		bool isMaxDepth = newDepth > 100;

		// Russian roulette for path termination
		bool isUseRR = newDepth > 5;
		bool isRR = isUseRR && rand() < obj->maxC;

		if (isMaxDepth || (isUseRR && !isRR))
			return obj->e;
		else {
			Vec f = (isUseRR && isRR) ? obj->cc : obj->c;
			Vec x = r.o + r.d * t;
			Vec n = (x - obj->p).norm();
			Vec nl = n % r.d < 0 ? n : n * -1;

			if (obj->refl == DIFF) { // Ideal DIFFUSE reflection
				double r1 = 2 * M_PI * rand();
				double r2 = rand();
				double r2s = sqrt(r2);

				Vec w = nl;
				Vec wo = w.x < -0.1 || w.x > 0.1 ? Vec::YAxis : Vec::XAxis;
				Vec u = (wo ^ w).norm();
				Vec v = w ^ u;

				Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();

				return obj->e + f * radiance(Ray(x, d), newDepth, rand);
			}
			else if (obj->refl == SPEC) // Ideal SPECULAR reflection
				return obj->e + f * radiance(Ray(x, r.d - n * (2 * (n % r.d))), newDepth, rand);
			else { // Ideal dielectric REFRACTION
				Ray reflRay(x, r.d - n * (2 * (n % r.d)));
				bool into = n % nl > 0;  // Ray from outside going in?
				double nc = 1;
				double nt = 1.5;
				double nnt = into ? nc / nt : nt / nc;
				double ddn = r.d % nl;
				double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);

				if (cos2t < 0)  // Total internal reflection
					return obj->e + f * radiance(reflRay, newDepth, rand);
				else
				{
					Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
					double a = nt - nc;
					double b = nt + nc;
					double R0 = (a * a) / (b * b);
					double c = 1 - (into ? -ddn : tdir % n);
					double Re = R0 + (1 - R0) * c * c * c * c * c;
					double Tr = 1 - Re;
					double P = .25 + .5 * Re;
					double RP = Re / P;
					double TP = Tr / (1 - P);

					Vec result;
					if (newDepth > 2) {
						// Russian roulette and splitting for selecting reflection and/or refraction
						if (rand() < P)
							result = radiance(reflRay, newDepth, rand) * RP;
						else
							result = radiance(Ray(x, tdir), newDepth, rand) * TP;
					}
					else
						result = radiance(reflRay, newDepth, rand) * Re + radiance(Ray(x, tdir), newDepth, rand) * Tr;

					return obj->e + f * result;
				}
			}
		}
	}
}

int main(int argc, char *argv[]) {
	clock_t start = clock(); 

	const int w = 256;
	const int h = 256;

	const int samps = argc == 2 ? atoi(argv[1]) / 4 : 25; // # samples

	const Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
	const Vec cx(w * .5135 / h);
	const Vec cy = (cx ^ cam.d).norm() * .5135;
	Vec* c = new Vec[w * h];

#pragma omp parallel for schedule(dynamic, 1)       // OpenMP
	// Loop over image rows
	for (int y = 0; y < h; y++) {
		fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
		RandomLCG rand(y);

		// Loop cols
		for (unsigned short x=0; x<w; x++) {
			// 2x2 subpixel rows
			for (int sy = 0; sy < 2; sy++) {
				const int i = (h - y - 1) * w + x;

				// 2x2 subpixel cols
				for (int sx = 0; sx < 2; sx++) {        
					Vec r = Vec::Zero;
					for (int s = 0; s < samps; s++) {
						double r1 = 2 * rand();
						double r2 = 2 * rand();
						double dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
						double dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);

						Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
							    cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;

						r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, rand) * (1.0 / samps);
					}
					c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
				}
			}
		}
	}

	printf("\n%f sec\n", (float)(clock() - start) / CLOCKS_PER_SEC); 

	FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i = 0; i < w * h; i++)
		fprintf(f,"%d %d %d\n", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
	fclose(f);
}
