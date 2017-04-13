using System;
using System.IO;

namespace Smallpt
{
    class RandomLCG
    {
        private uint mSeed;

        public RandomLCG(uint seed)
        {
            mSeed = seed;
        }

        public double NextDouble()
        {
            mSeed = 214013u * mSeed + 2531011u;
            return mSeed * (1.0 / 4294967296.0);
        }
    }

    struct Vec
    {
        public double x;
        public double y;
        public double z;

        public static Vec Zero = new Vec(0, 0, 0);
        public static Vec XAxis = new Vec(1, 0, 0);
        public static Vec YAxis = new Vec(0, 1, 0);
        public static Vec ZAxis = new Vec(0, 0, 1);

        public Vec(double x, double y, double z)
        {
            this.x = x; 
            this.y = y;
            this.z = z;
        }

        public static Vec Add(ref Vec a, ref Vec b)
        {
            return new Vec(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        public static Vec Sub(ref Vec a, ref Vec b)
        {
            return new Vec(a.x - b.x, a.y - b.y, a.z - b.z);
        }

        public static Vec Mul(ref Vec a, double b)
        {
            return new Vec(a.x * b, a.y * b, a.z * b);
        }

        // component-wise multiplication
        public static Vec mult(ref Vec a, ref Vec b)
        {
            return new Vec(a.x * b.x, a.y * b.y, a.z * b.z);
        }

        public Vec norm()
        {
            double s = 1.0 / Math.Sqrt(x * x + y * y + z * z);
            return new Vec(x * s, y * s, z * s);
        }

        public static double dot(ref Vec a, ref Vec b)
        {
            return a.x * b.x + a.y * b.y + a.z * b.z;
        }

        // cross product
        public static Vec cross(ref Vec a, ref Vec b)
        {
            return new Vec(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
        }

        public static Vec operator +(Vec a, Vec b)
        {
            return new Vec(a.x + b.x, a.y + b.y, a.z + b.z);
        }

        public static Vec operator -(Vec a, Vec b)
        {
            return new Vec(a.x - b.x, a.y - b.y, a.z - b.z);
        }

        public static Vec operator *(Vec a, double b)
        {
            return new Vec(a.x * b, a.y * b, a.z * b);
        }

        // component-wise multiplication
        public Vec mult(Vec b)
        {
            return new Vec(x * b.x, y * b.y, z * b.z);
        }

        // cross product
        public static Vec operator %(Vec a, Vec b)
        {
            return new Vec(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
        }
    }

    // material types
    enum Refl_t
    {
        DIFF,
        SPEC,
        REFR
    };

    class Ray
    {
        public Vec o;
        public Vec d;

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
            maxC = Math.Max(Math.Max(c.x, c.y), c.z);
            cc = Vec.Mul(ref c, 1.0 / maxC);
        }

        // returns distance, 1e20 if nohit
        public double intersect(Ray r)
        {
            // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
            Vec op = Vec.Sub(ref p, ref r.o);
            double b = Vec.dot(ref op, ref r.d);
            double det = b * b - Vec.dot(ref op, ref op) + sqRad;
            const double eps = 1e-4;

            if (det < 0)
                return 1e20;
            else
            {
                double dets = Math.Sqrt(det);

                if (b - dets > eps)
                    return b - dets;
                else if (b + dets > eps)
                    return b + dets;
                else
                    return 1e20;
            }
        }
    };

    class Smallpt
    {
        //Scene: radius, position, emission, color, material
        static Vec a = new Vec(1,1,1);
        static Sphere[] spheres =
        {
			new Sphere(1e5,  new Vec( 1e5+1,40.8,81.6),  Vec.Zero, new Vec(.75,.25,.25), Refl_t.DIFF),//Left
			new Sphere(1e5,  new Vec(-1e5+99,40.8,81.6), Vec.Zero, new Vec(.25,.25,.75), Refl_t.DIFF),//Rght
			new Sphere(1e5,  new Vec(50,40.8, 1e5),      Vec.Zero, new Vec(.75,.75,.75), Refl_t.DIFF),//Back
			new Sphere(1e5,  new Vec(50,40.8,-1e5+170),  Vec.Zero, Vec.Zero,             Refl_t.DIFF),//Frnt
			new Sphere(1e5,  new Vec(50, 1e5, 81.6),     Vec.Zero, new Vec(.75,.75,.75), Refl_t.DIFF),//Botm
			new Sphere(1e5,  new Vec(50,-1e5+81.6,81.6), Vec.Zero, new Vec(.75,.75,.75), Refl_t.DIFF),//Top
			new Sphere(16.5, new Vec(27,16.5,47),        Vec.Zero, Vec.Mul(ref a, .999),  Refl_t.SPEC),//Mirr
			new Sphere(16.5, new Vec(73,16.5,78),        Vec.Zero, Vec.Mul(ref a, .999),  Refl_t.REFR),//Glas
			new Sphere(600,  new Vec(50,681.6-.27,81.6), new Vec(12,12,12), Vec.Zero,    Refl_t.DIFF) //Lite
		};

        //static Random random = new Random();
        static RandomLCG random = new RandomLCG(0u);

        static double rand()
        {
            return random.NextDouble();
        }

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
            return (int)(Math.Pow(clamp(x), 1 / 2.2) * 255 + .5);
        }

        static Sphere intersect(Ray r, out double t)
        {
            double d, inf = t = 1e20;
            Sphere ret = null;

            foreach (Sphere s in spheres)
            {
                d = s.intersect(r);
                if (d < t)
                {
                    t = d;
                    ret = s;
                }
            }

            return ret;
        }

        static Vec radiance(Ray r, int depth)
        {
            double t;   // distance to intersection
            Sphere obj = intersect(r, out t);

            if (obj == null)
                return Vec.Zero;       // if miss, return black
            else
            {
                int newDepth = depth + 1;
                bool isMaxDepth = newDepth > 100;

                // Russian roulette for path termination
                bool isUseRR = newDepth > 5;
                bool isRR = isUseRR && rand() < obj.maxC;

                if (isMaxDepth || (isUseRR && !isRR))
                    return obj.e;
                else
                {
                    Vec f = (isUseRR && isRR) ? obj.cc : obj.c;
                    Vec aa = Vec.Mul(ref r.d, t);
                    Vec x = Vec.Add(ref r.o, ref aa);
                    Vec n = Vec.Sub(ref x, ref obj.p).norm();
                    Vec nl = Vec.dot(ref n, ref r.d) < 0 ? n : Vec.Mul(ref n, -1);

                    if (obj.refl == Refl_t.DIFF) // Ideal DIFFUSE reflection
                    {
                        double r1 = 2 * Math.PI * rand();
                        double r2 = rand();
                        double r2s = Math.Sqrt(r2);

                        Vec w = nl;
                        Vec wo = Math.Abs(w.x) > .1 ? Vec.YAxis : Vec.XAxis;
                        Vec u = Vec.cross(ref wo, ref w).norm();
                        Vec v = Vec.cross(ref w, ref u);

                        Vec d1 = Vec.Mul(ref u, Math.Cos(r1) * r2s);
                        Vec d2 = Vec.Mul(ref v, Math.Sin(r1) * r2s);
                        Vec d3 = Vec.Mul(ref w,  Math.Sqrt(1 - r2));
                        Vec d4 = Vec.Add(ref d1, ref d2);
                        Vec d =  Vec.Add(ref d4, ref d3).norm();

                        Vec k1 = radiance(new Ray(x, d), newDepth);
                        Vec k2 = Vec.mult(ref f, ref k1);
                        return Vec.Add(ref obj.e, ref k2);
                    }
                    else if (obj.refl == Refl_t.SPEC) // Ideal SPECULAR reflection
                        return obj.e + f.mult(radiance(new Ray(x, r.d - n * 2 * Vec.dot(ref n, ref r.d)), newDepth));
                    else // Ideal dielectric REFRACTION
                    {
                        Ray reflRay = new Ray(x, r.d - n * 2 * Vec.dot(ref n, ref r.d));
                        bool into = Vec.dot(ref n, ref nl) > 0;  // Ray from outside going in?
                        double nc = 1;
                        double nt = 1.5;
                        double nnt = into ? nc / nt : nt / nc;
                        double ddn = Vec.dot(ref r.d, ref nl);
                        double cos2t = 1 - nnt * nnt * (1 - ddn * ddn);

                        if (cos2t < 0)  // Total internal reflection
                            return obj.e + f.mult(radiance(reflRay, newDepth));
                        else
                        {
                            Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + Math.Sqrt(cos2t)))).norm();
                            double a = nt - nc;
                            double b = nt + nc;
                            double R0 = (a * a) / (b * b);
                            double c = 1 - (into ? -ddn : Vec.dot(ref tdir, ref n));
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
                                    result = radiance(reflRay, newDepth) * RP;
                                else
                                    result = radiance(new Ray(x, tdir), newDepth) * TP;
                            }
                            else
                                result = radiance(reflRay, newDepth) * Re + radiance(new Ray(x, tdir), newDepth) * Tr;

                            return obj.e + f.mult(result);
                        }
                    }
                }
            }
        }

        public static void Main(string[] args)
        {
            DateTime start = DateTime.Now;

            const int w = 256;
            const int h = 256;

            int samps = args.Length == 2 ? int.Parse(args[1]) / 4 : 25; // # samples

            // cam pos, dir
            Ray cam = new Ray(new Vec(50, 52, 295.6), new Vec(0, -0.042612, -1).norm());
            Vec cx = new Vec(w * .5135 / h, 0, 0);
            Vec cy = (cx % cam.d).norm() * .5135;

            // final color buffer
            Vec[] c = new Vec[w * h];

            // Loop over image rows
            for (int y = 0; y < h; y++)
            {
                Console.Write("\rRendering ({0} spp) {1:F2}%", samps * 4, 100.0 * y / (h - 1));

                // Loop cols
                for (int x = 0; x < w; x++)
                {
                    int i = (h - y - 1) * w + x;
                    c[i] = Vec.Zero;

                    // 2x2 subpixel rows
                    for (int sy = 0; sy < 2; sy++)
                    {
                        // 2x2 subpixel cols
                        for (int sx = 0; sx < 2; sx++)
                        {
                            Vec r = Vec.Zero;
                            for (int s = 0; s < samps; s++)
                            {
                                double r1 = 2 * rand();
                                double r2 = 2 * rand();
                                double dx = r1 < 1 ? Math.Sqrt(r1) - 1 : 1 - Math.Sqrt(2 - r1);
                                double dy = r2 < 1 ? Math.Sqrt(r2) - 1 : 1 - Math.Sqrt(2 - r2);
                                Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                        cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;

                                // Camera rays are pushed forward to start in interior
                                Ray camRay = new Ray(cam.o + d * 140, d.norm());

                                // Accumuate radiance
                                r = r + radiance(camRay, 0) * (1.0 / samps);
                            }

                            // Convert radiance to color
                            c[i] = c[i] + new Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                        }
                    }
                }
            }

            Console.WriteLine("\n{0} sec", (DateTime.Now - start).TotalSeconds);

            using (StreamWriter sw = new StreamWriter("image.ppm"))
            {
                sw.Write("P3\r\n{0} {1}\r\n{2}\r\n", w, h, 255);
                for (int i = 0; i < w * h; i++)
                    sw.Write("{0} {1} {2}\r\n", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
            }
        }
    }
}