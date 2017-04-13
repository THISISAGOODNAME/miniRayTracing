import java.io.*;

class RandomLCG
{
	long mSeed; // java has no unsigned int
	
	public RandomLCG(long seed) 
	{
		mSeed = seed;
	}
	
	public double nextDouble()
	{ 
		mSeed = (214013L * mSeed + 2531011L) % 0x100000000L;
		return mSeed * (1.0 / 4294967296.0);
	}
}

class Vec
{
	public final double x, y, z;	// position, also color (r,g,b)
	public static final Vec Zero = new Vec(0, 0, 0);
	
	public Vec(double x, double y, double z) 
	{
		this.x = x; 
		this.y = y; 
		this.z = z; 
	}

	public Vec add(Vec b) 
	{
		return new Vec(x + b.x, y + b.y, z + b.z); 
	}

	public Vec sub(Vec b) 
	{
		return new Vec(x - b.x, y - b.y, z - b.z); 
	}

	public Vec mul(double b) 
	{
		return new Vec(x * b, y * b, z * b); 
	}

	public Vec mult(Vec b) 
	{
		return new Vec(x * b.x, y * b.y, z * b.z); 
	}

	public Vec norm() 
	{
		return this.mul(1 / Math.sqrt(x * x + y * y + z * z)); 
	}

	public double dot(Vec b) 
	{
		return x * b.x + y * b.y + z * b.z; 
	}

	public Vec cross(Vec b) 
	{
		return new Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); 
	}
}

// material types, used in radiance()
enum Refl_t
{ 
	DIFF, 
	SPEC, 
	REFR 
};

class Ray
{
	public Vec o, d;
	
	public Ray(Vec o, Vec d)
	{ 
		this.o = o; 
		this.d = d; 
	}
}

class Sphere 
{
	public double rad;       // radius
	public Vec p, e, c;      // position, emission, color
	public Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
	public double maxC;
	public Vec cc;
	private double sqRad;
	
	public Sphere(double rad, Vec p, Vec e, Vec c, Refl_t refl)
	{
		this.rad = rad; 
		this.p = p; 
		this.e = e; 
		this.c = c; 
		this.refl = refl; 
		
		sqRad = rad * rad;
		maxC = Math.max(Math.max(c.x, c.y), c.z);
		cc = c.mul(1.0 / maxC);
	}
	
	// returns distance, 0 if nohit
	public double intersect(Ray r) 
	{ 
		Vec op = p.sub(r.o); // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double eps = 1e-4;
		double b = op.dot(r.d);
		double det = b * b - op.dot(op) + sqRad;
		
		if (det < 0) 
			return 0; 
		else 
		{
			double dets = Math.sqrt(det);

			if (b - dets > eps)
				return b - dets;
			else if (b + dets > eps)
				return b + dets;
			else
				return 0;
		}
	}
}

// Java has no pass by ref
class IntersectResult
{ 
	public double t = 1e20;  
	public Sphere s;
}

class smallpt 
{
	//static java.util.Random random = new java.util.Random();
	static RandomLCG random = new RandomLCG(0L);
	
	static double rand()
	{ 
		return random.nextDouble(); 
	}
	
	//Scene: radius, position, emission, color, material
	static Sphere[] spheres = 
	{
		new Sphere(1e5, new Vec( 1e5+1,40.8,81.6), Vec.Zero,new Vec(.75,.25,.25),Refl_t.DIFF),//Left
		new Sphere(1e5, new Vec(-1e5+99,40.8,81.6),Vec.Zero,new Vec(.25,.25,.75),Refl_t.DIFF),//Rght
		new Sphere(1e5, new Vec(50,40.8, 1e5),     Vec.Zero,new Vec(.75,.75,.75),Refl_t.DIFF),//Back
		new Sphere(1e5, new Vec(50,40.8,-1e5+170), Vec.Zero,Vec.Zero,           Refl_t.DIFF),//Frnt
		new Sphere(1e5, new Vec(50, 1e5, 81.6),    Vec.Zero,new Vec(.75,.75,.75),Refl_t.DIFF),//Botm
		new Sphere(1e5, new Vec(50,-1e5+81.6,81.6),Vec.Zero,new Vec(.75,.75,.75),Refl_t.DIFF),//Top
		new Sphere(16.5,new Vec(27,16.5,47),       Vec.Zero,new Vec(1,1,1).mul(.999), Refl_t.SPEC),//Mirr
		new Sphere(16.5,new Vec(73,16.5,78),       Vec.Zero,new Vec(1,1,1).mul(.999), Refl_t.REFR),//Glas
		new Sphere(600, new Vec(50,681.6-.27,81.6),new Vec(12,12,12),  Vec.Zero, Refl_t.DIFF) //Lite
	};
	
	static double clamp(double x)
	{ 
		if (x < 0)
			return 0;
		else if (x > 1)
			return 1;
		else
			return x;
	}
	
	static int toInt(double x) 
	{ 
		return (int)(Math.pow(clamp(x), 1 / 2.2) * 255 + .5); 
	}
	
	static IntersectResult intersect(Ray r) 
	{
		IntersectResult result = new IntersectResult();
		
		for (Sphere s : spheres)
		{
			double d = s.intersect(r);
			if (d != 0 && d < result.t) 
			{ 
				result.t = d; 
				result.s = s; 
			}
		}
		return result;
	}
	
	static Vec radiance(Ray r, int depth) 
	{
		IntersectResult ires = intersect(r);
	   
		if (ires.s == null) 
		{
			return Vec.Zero; // if miss, return black
		}
		else
		{
			Sphere obj = ires.s;        // the hit object
			double t = ires.t;
			
			int newDepth = depth + 1;
			boolean isMaxDepth = newDepth > 100;

			// Russian roulette for path termination
			boolean isUseRR = newDepth > 5;
			boolean isRR = isUseRR && rand() < obj.maxC;

			if (isMaxDepth || (isUseRR && !isRR))
				return obj.e;
			else
			{
				Vec f = (isUseRR && isRR) ? obj.cc : obj.c;
				Vec x = r.o.add(r.d.mul(t));
				Vec n = x.sub(obj.p).norm();
				Vec nl = n.dot(r.d) < 0 ? n : n.mul(-1);
			
				if (obj.refl == Refl_t.DIFF) // Ideal DIFFUSE reflection
				{
					double r1 = 2 * Math.PI * rand();
					double r2 = rand();
					double r2s = Math.sqrt(r2);
					
					Vec w = nl;
					Vec wo = Math.abs(w.x) > .1 ? new Vec(0, 1, 0) : new Vec(1, 0, 0);
					Vec u = wo.cross(w).norm();
					Vec v = w.cross(u);

					Vec d = (u.mul(Math.cos(r1) * r2s).add(v.mul(Math.sin(r1) * r2s)).add(w.mul(Math.sqrt(1 - r2)))).norm();

					return obj.e.add(f.mult(radiance(new Ray(x, d), newDepth)));
				}
				else if (obj.refl == Refl_t.SPEC) // Ideal SPECULAR reflection
				{
					return obj.e.add(f.mult(radiance(new Ray(x, r.d.sub(n.mul(2 * n.dot(r.d)))), newDepth)));
				}
				else // Ideal dielectric REFRACTION
				{
					Ray reflRay = new Ray(x, r.d.sub(n.mul(2 * n.dot(r.d))));     
					boolean into = n.dot(nl) > 0;	// Ray from outside going in?
					double nc = 1;
					double nt = 1.5;
					double nnt = into ? nc / nt : nt / nc;
					double ddn = r.d.dot(nl);
					double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);
					
					if (cos2t < 0)	// Total internal reflection
					{
						return obj.e.add(f.mult(radiance(reflRay, newDepth)));
					}
					else
					{
						Vec tdir = (r.d.mul(nnt).sub(n.mul((into ? 1 : -1) * (ddn * nnt + Math.sqrt(cos2t))))).norm();
						double a = nt - nc;
						double b = nt + nc;
						double R0 = a * a / (b * b);
						double c = 1 - (into ? -ddn : tdir.dot(n));
						double Re = R0 + (1 - R0) * c * c * c * c * c;
						double Tr = 1 - Re;
						double P = .25 + .5 * Re;
						double RP = Re / P;
						double TP = Tr / (1 - P);
					
						Vec result;
						if (newDepth > 2)
						{
							// Russian roulette and splitting for selecting reflection and/or refraction
							if (rand() < P)
								result = radiance(reflRay, newDepth).mul(RP);
							else
								result = radiance(new Ray(x, tdir), newDepth).mul(TP);
						}
						else
							result = radiance(reflRay, newDepth).mul(Re).add(radiance(new Ray(x, tdir), newDepth).mul(Tr));

						return obj.e.add(f.mult(result));
					}
				}
			}
		}
	}
	
	public static void main(String[] args)
	{
		long start = System.currentTimeMillis();
		
		int w = 256;
		int h = 256;
		int samps = args.length == 2 ? Integer.parseInt(args[1]) / 4 : 25; // # samples
		
		// cam pos, dir
		Ray cam = new Ray(new Vec(50, 52, 295.6), new Vec(0, -0.042612, -1).norm());
		Vec cx = new Vec(w * .5135 / h, 0, 0);
		Vec cy = (cx.cross(cam.d)).norm().mul(.5135);
		
		// final color buffer
		Vec[] c = new Vec[w * h];
			
		// Loop over image rows
		for (int y = 0; y < h; y++)
		{
			System.out.printf("\rRendering (%d spp) %5.2f%%", samps * 4, 100.0 * y / (h - 1));
			
			// Loop cols
			for (int x = 0; x < w; x++)
			{
				// 2x2 subpixel rows
				for (int sy = 0; sy < 2; sy++)
				{
					int i = (h - y - 1) * w + x;
					c[i] = Vec.Zero;
					
					// 2x2 subpixel cols
					for (int sx = 0; sx < 2; sx++)
					{
						Vec r = Vec.Zero;
						for (int s = 0; s < samps; s++)
						{
							double r1 = 2 * rand();
							double r2 = 2 * rand();
							double dx = r1 < 1 ? Math.sqrt(r1) - 1 : 1 - Math.sqrt(2 - r1);
							double dy = r2 < 1 ? Math.sqrt(r2) - 1 : 1 - Math.sqrt(2 - r2);

							Vec d = cx.mul(((sx + .5 + dx) / 2 + x) / w - .5).add(
									cy.mul(((sy + .5 + dy) / 2 + y) / h - .5)).add(cam.d);

							// Camera rays are pushed forward to start in interior
							Ray camRay = new Ray(cam.o.add(d.mul(140)), d.norm());

							// Accumuate radiance
							r = r.add(radiance(camRay, 0).mul(1.0 / samps));
						}
						
						// Convert radiance to color
						c[i] = c[i].add((new Vec(clamp(r.x), clamp(r.y), clamp(r.z))).mul(.25));
					}
				}
			}
		}
		
		System.out.printf("\n%f sec", (System.currentTimeMillis() - start) / 1000.0);
		
		try
		{
			FileWriter fw = new FileWriter("image.ppm");
			fw.write("P3\r\n" + w + " " + h + "\r\n255\r\n");
			for (int i = 0; i < w * h; i++)
				fw.write(toInt(c[i].x) + " " + toInt(c[i].y) + " " + toInt(c[i].z) + "\r\n");
			fw.close();
		}
		catch(java.io.IOException e) 
		{
		}
	}
}
