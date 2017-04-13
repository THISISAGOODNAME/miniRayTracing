import math
import time
import sys

class RandomLCG:
	def __init__(self, seed):
		self.seed = seed

	def __call__(self) :
		self.seed = (214013 * self.seed + 2531011) % 0x100000000
		return self.seed * (1.0 / 4294967296.0)

class Vec:
	def __init__(self, x, y, z):
		self.x = x
		self.y = y
		self.z = z

	def __add__(self, b):
		return Vec(self.x + b.x, self.y + b.y, self.z + b.z)

	def __sub__(self, b):
		return Vec(self.x - b.x, self.y - b.y, self.z - b.z)

	def __mul__(self, b):
		return Vec(self.x * b, self.y * b, self.z * b)

	def mult(self, b):
		return Vec(self.x * b.x, self.y * b.y, self.z * b.z)

	def norm(self):
		return self * (1.0 / math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z))

	def __pow__(self, b): # dot product
		return self.x * b.x + self.y * b.y + self.z * b.z

	def __mod__(self, b):
		return Vec(self.y * b.z - self.z * b.y, self.z * b.x - self.x * b.z, self.x * b.y - self.y * b.x)

Vec.Zero = Vec(0, 0, 0)
Vec.XAxis = Vec(1, 0, 0)
Vec.YAxis = Vec(0, 1, 0)
Vec.ZAxis = Vec(0, 0, 1)

class Refl:
	(DIFF, SPEC, REFR) = range(3)

class Ray:
	def __init__(self, o, d):
		self.o = o
		self.d = d

class Sphere:
	def __init__(self, rad, p, e, c, refl):
		self.rad = rad
		self.p = p
		self.e = e
		self.c = c
		self.refl = refl
		self.sqRad = rad * rad
		self.maxC = max(c.x, c.y, c.z)
		if self.maxC == 0:
			self.cc = Vec.Zero
		else:
			self.cc = c * (1.0 / self.maxC)

	def intersect(self, r):
		# Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		op = self.p - r.o
		b = op ** r.d
		det = b * b - op ** op + self.sqRad
		eps = 1e-4

		if det < 0:
			return 0
		else:
			dets = math.sqrt(det)

			if b - dets > eps:
				return b - dets
			elif b + dets > eps:
				return b + dets
			else:
				return 0

# Scene: radius, position, emission, color, material
spheres = [
	Sphere(1e5,  Vec( 1e5+1,40.8,81.6),  Vec.Zero, Vec(.75,.25,.25), Refl.DIFF), #Left
	Sphere(1e5,  Vec(-1e5+99,40.8,81.6), Vec.Zero, Vec(.25,.25,.75), Refl.DIFF), #Rght
	Sphere(1e5,  Vec(50,40.8, 1e5),      Vec.Zero, Vec(.75,.75,.75), Refl.DIFF), #Back
	Sphere(1e5,  Vec(50,40.8,-1e5+170),  Vec.Zero, Vec.Zero,         Refl.DIFF), #Frnt
	Sphere(1e5,  Vec(50, 1e5, 81.6),     Vec.Zero, Vec(.75,.75,.75), Refl.DIFF), #Botm
	Sphere(1e5,  Vec(50,-1e5+81.6,81.6), Vec.Zero, Vec(.75,.75,.75), Refl.DIFF), #Top
	Sphere(16.5, Vec(27,16.5,47),        Vec.Zero, Vec(1,1,1)*.999,  Refl.SPEC), #Mirr
	Sphere(16.5, Vec(73,16.5,78),        Vec.Zero, Vec(1,1,1)*.999,  Refl.REFR), #Glas
	Sphere(600,  Vec(50,681.6-.27,81.6), Vec(12,12,12), Vec.Zero,    Refl.DIFF)  #Lite
]

rand = RandomLCG(0)

def clamp(x):
	if x < 0:
		return 0 
	elif x > 1:
		return 1
	else:
		return x 

def toInt(x):
	return int((clamp(x) ** (1 / 2.2)) * 255 + .5)

def intersect(r):
	t = 1e20
	for s in spheres:
		d = s.intersect(r)
		if d != 0 and d < t:
			t = d
			obj = s
	return obj, t

def radiance(r, depth):
	obj, t = intersect(r)

	if obj == None: 
		return Vec.Zero       # if miss, return black
	else:
		newDepth = depth + 1
		isMaxDepth = newDepth > 100

		# Russian roulette for path termination
		isUseRR = newDepth > 5
		isRR = isUseRR and rand() < obj.maxC

		if isMaxDepth or (isUseRR and not isRR):
			return obj.e
		else:
			f = (isUseRR and isRR) and obj.cc or obj.c
			x = r.o + r.d * t
			n = (x - obj.p).norm()
			nl = n ** r.d < 0 and n or n * -1

			if obj.refl == Refl.DIFF: # Ideal DIFFUSE reflection
				r1 = 2 * math.pi * rand()
				r2 = rand()
				r2s = math.sqrt(r2)

				w = nl
				wo = abs(w.x) > .1 and Vec.YAxis or Vec.XAxis
				u = (wo % w).norm()
				v = w % u

				d = (u * (math.cos(r1) * r2s) + v * (math.sin(r1) * r2s) + w * math.sqrt(1 - r2)).norm()

				return obj.e + f.mult(radiance(Ray(x, d), newDepth))
			
			elif obj.refl == Refl.SPEC: # Ideal SPECULAR reflection
				return obj.e + f.mult(radiance(Ray(x, r.d - n * (2 * n ** r.d)), newDepth))
			
			else: # Ideal dielectric REFRACTION
				reflRay = Ray(x, r.d - n * (2 * n ** r.d))
				into = n ** nl > 0 # from outside going in?
				nc = 1
				nt = 1.5
				nnt = into and nc / nt or nt / nc
				ddn = r.d ** nl
				cos2t = 1 - nnt * nnt * (1 - ddn * ddn)

				if cos2t < 0:	# Total internal reflection
				
					return obj.e + f.mult(radiance(reflRay, newDepth))
				
				else:
					tdir = (r.d * nnt - n * ((into and 1 or -1) * (ddn * nnt + math.sqrt(cos2t)))).norm()
					a = nt - nc
					b = nt + nc
					R0 = a * a / (b * b)
					c = 1 - (into and  -ddn or tdir ** n)
					Re = R0 + (1 - R0) * c * c * c * c * c
					Tr = 1 - Re
					P = .25 + .5 * Re
					RP = Re / P
					TP = Tr / (1 - P)

					if newDepth > 2:
						# Russian roulette and splitting for selecting reflection and/or refraction
						if rand() < P:
							result = radiance(reflRay, newDepth) * RP
						else:
							result = radiance(Ray(x, tdir), newDepth) * TP
					else:
						result = radiance(reflRay, newDepth)* Re + radiance(Ray(x, tdir), newDepth) * Tr

					return obj.e + f.mult(result)

start = time.clock()

w = 256
h = 256

if len(sys.argv) == 2:
	samps = int(sys.argv[1]) / 4
else:
	samps = 25

# cam pos, dir
cam = Ray(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm())
cx = Vec(w * .5135 / h, 0, 0)
cy = (cx % cam.d).norm() * .5135

# final color buffer
c = {}

# Loop over image rows
for y in range(0, h):
	sys.stderr.write('\rRendering ({0} spp) {1:2.2%}'.format(samps * 4, y / (h - 1)))

	# Loop cols
	for x in range(0, w):
		i = (h - y - 1) * w + x
		c[i] = Vec.Zero

		# 2x2 subpixel rows
		for sy in range (0, 2):
			# 2x2 subpixel cols
			for sx in range (0, 2):
				r = Vec.Zero
				for s in range(samps):
					r1 = 2 * rand()
					r2 = 2 * rand()
					dx = (r1 < 1) and (math.sqrt(r1) - 1) or (1 - math.sqrt(2 - r1))
					dy = (r2 < 1) and (math.sqrt(r2) - 1) or (1 - math.sqrt(2 - r2))
					d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) + \
						cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d

					# Camera rays are pushed forward to start in interior
					camRay = Ray(cam.o + d * 140, d.norm())

					# Accumuate radiance
					r = r + radiance(camRay, 0) * (1.0 / samps)

				# Convert radiance to color
				c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25

print('\n{0} sec'.format(time.clock() - start))

f = open('image.ppm', 'w')
f.write('P3\n{0} {1}\n{2}\n'.format(w, h, 255))
for i in range(0, w * h):
	f.write('{0} {1} {2}\n'.format(toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)))

