class RandomLCG
	def initialize(seed)
		@seed = seed
	end

	def next
		@seed = (214013 * @seed + 2531011) % 0x100000000
		@seed * (1.0 / 4294967296.0)
	end
end

class Vec
	attr_reader :x, :y, :z
	
	def initialize(x, y, z)
		@x = x
		@y = y
		@z = z
	end

	def +(b)
		Vec.new(@x + b.x, @y + b.y, @z + b.z)
	end

	def -(b)
		Vec.new(@x - b.x, @y - b.y, @z - b.z)
	end

	def *(b)
		Vec.new(@x * b, @y * b, @z * b)
	end

	def mult(b)
		Vec.new(@x * b.x, @y * b.y, @z * b.z)
	end

	def norm
		self * (1.0 / Math.sqrt(@x * @x + @y * @y + @z * @z))
	end

	# dot product
	def **(b)
		@x * b.x + @y * b.y + @z * b.z
	end

	# cross product
	def %(b)
		Vec.new(@y * b.z - @z * b.y, @z * b.x - @x * b.z, @x * b.y - @y * b.x)
	end

	Zero = Vec.new(0.0, 0.0, 0.0)
	XAxis = Vec.new(1.0, 0.0, 0.0)
	YAxis = Vec.new(0.0, 1.0, 0.0)
	ZAxis = Vec.new(0.0, 0.0, 1.0)
end

class Refl
	DIFF = 1
	SPEC = 2
	REFR = 3
end

class Ray
	attr_reader :o, :d
	def initialize(o, d)
		@o = o
		@d = d
	end
end

class Sphere
	attr_reader :rad, :p, :e, :c, :refl, :maxC, :cc

	def initialize(rad, p, e, c, refl)
		@rad = rad
		@p = p
		@e = e
		@c = c
		@refl = refl
		@sqRad = rad * rad
		@maxC = [c.x, c.y, c.z].max
		@cc = @maxC == 0 ? Vec::Zero : c * (1.0 / @maxC)
	end

	def intersect(r)
		# Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		op = @p - r.o
		b = op ** r.d
		det = b * b - op ** op + @sqRad
		eps = 1e-4

		if det < 0
			0
		else
			dets = Math.sqrt(det)

			if b - dets > eps
				b - dets
			elsif b + dets > eps
				b + dets
			else
				0
			end
		end
	end
end

$spheres = [
	Sphere.new(1e5,   Vec.new( 1e5+1,40.8,81.6),   Vec::Zero, Vec.new(0.75,0.25,0.25),     Refl::DIFF), #Left
	Sphere.new(1e5,   Vec.new(-1e5+99,40.8,81.6),  Vec::Zero, Vec.new(0.25,0.25,0.75),     Refl::DIFF), #Rght
	Sphere.new(1e5,   Vec.new(50.0,40.8, 1e5),      Vec::Zero, Vec.new(0.75,0.75,0.75),    Refl::DIFF), #Back
	Sphere.new(1e5,   Vec.new(50.0,40.8,-1e5+170),  Vec::Zero, Vec::Zero,				   Refl::DIFF), #Frnt
	Sphere.new(1e5,   Vec.new(50.0, 1e5, 81.6),     Vec::Zero, Vec.new(0.75,0.75,0.75),    Refl::DIFF), #Botm
	Sphere.new(1e5,   Vec.new(50.0,-1e5+81.6,81.6), Vec::Zero, Vec.new(0.75,0.75,0.75),    Refl::DIFF), #Top
	Sphere.new(16.5,  Vec.new(27.0,16.5,47),        Vec::Zero, Vec.new(1.0,1.0,1.0)*0.999, Refl::SPEC), #Mirr
	Sphere.new(16.5,  Vec.new(73.0,16.5,78),        Vec::Zero, Vec.new(1.0,1.0,1.0)*0.999, Refl::REFR), #Glas
	Sphere.new(600.0, Vec.new(50,681.6-0.27,81.6),  Vec.new(12.0,12.0,12.0), Vec::Zero,    Refl::DIFF)  #Lite
]

$random = RandomLCG.new(0)
def rand
	$random.next
end

def clamp(x)
	if x < 0.0
		0.0
	elsif x > 1
		1.0
	else
		x
	end
end

def toInt(x)
	((clamp(x) ** (1.0 / 2.2)) * 255.0 + 0.5).truncate
end

def intersect(r)
	t = 1e20
	for s in $spheres
		d = s.intersect(r)
		if d != 0 and d < t
			t = d
			obj = s
		end
	end
	return obj, t
end

def radiance(r, depth)
	obj, t = intersect(r)

	if obj == nil
		return Vec::Zero       # if miss, return black
	else
		newDepth = depth + 1
		isMaxDepth = newDepth > 100

		# Russian roulette for path termination
		isUseRR = newDepth > 5
		isRR = isUseRR and rand < obj.maxC

		if isMaxDepth or (isUseRR and not isRR)
			return obj.e
		else
			f = (isUseRR and isRR) ? obj.cc : obj.c
			x = r.o + r.d * t
			n = (x - obj.p).norm
			nl = n ** r.d < 0 ? n : n * -1.0

			if obj.refl == Refl::DIFF # Ideal DIFFUSE reflection
				r1 = 2 * Math::PI * rand
				r2 = rand
				r2s = Math.sqrt(r2)

				w = nl
				wo = w.x.abs > 0.1 ? Vec::YAxis : Vec::XAxis
				u = (wo % w).norm
				v = w % u

				d = (u * (Math.cos(r1) * r2s) + v * (Math.sin(r1) * r2s) + w * Math.sqrt(1 - r2)).norm

				obj.e + f.mult(radiance(Ray.new(x, d), newDepth))
			
			elsif obj.refl == Refl::SPEC # Ideal SPECULAR reflection
				obj.e + f.mult(radiance(Ray.new(x, r.d - n * (2 * n ** r.d)), newDepth))
			
			else # Ideal dielectric REFRACTION
				reflRay = Ray.new(x, r.d - n * (2.0 * n ** r.d))
				into = n ** nl > 0 # from outside going in?
				nc = 1.0
				nt = 1.5
				nnt = into ? nc / nt : nt / nc
				ddn = r.d ** nl
				cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn)

				if cos2t < 0	# Total internal reflection
					obj.e + f.mult(radiance(reflRay, newDepth))
				else
					tdir = (r.d * nnt - n * ((into and 1 or -1) * (ddn * nnt + Math.sqrt(cos2t)))).norm
					a = nt - nc
					b = nt + nc
					r0 = a * a / (b * b)
					c = 1.0 - (into and  -ddn or tdir ** n)
					re = r0 + (1 - r0) * c * c * c * c * c
					tr = 1.0 - re
					p = 0.25 + 0.5 * re
					rp = re / p
					tp = tr / (1 - p)

					if newDepth > 2
						# Russian roulette and splitting for selecting reflection and/or refraction
						if rand < p
							result = radiance(reflRay, newDepth) * rp
						else
							result = radiance(Ray.new(x, tdir), newDepth) * tp
						end
					else
						result = radiance(reflRay, newDepth)* re + radiance(Ray.new(x, tdir), newDepth) * tr
					end

					obj.e + f.mult(result)
				end
			end
		end
	end
end

start = Time.now

w = 256
h = 256

if ARGV.length == 1
	samps = ARGV[0].to_i / 4
else
	samps = 1
end

# cam pos, dir
cam = Ray.new(Vec.new(50.0, 52.0, 295.6), Vec.new(0.0, -0.042612, -1.0).norm)
cx = Vec.new(w * 0.5135 / h, 0.0, 0.0)
cy = (cx % cam.d).norm * 0.5135

# final color buffer
c = {}

# Loop over image rows
for y in 0..h-1
	$stderr.print "\rRendering (#{samps * 4} spp) #{"%2.2f" % (100.0 * y / (h - 1))}%"

	# Loop cols
	for x in 0..w-1
		i = (h - y - 1) * w + x
		c[i] = Vec::Zero

		# 2x2 subpixel rows
		for sy in 0..1
			# 2x2 subpixel cols
			for sx in 0..1
				r = Vec::Zero
				for s in 1..samps
					r1 = 2 * rand
					r2 = 2 * rand
					dx = (r1 < 1) ? (Math.sqrt(r1) - 1) : (1 - Math.sqrt(2 - r1))
					dy = (r2 < 1) ? (Math.sqrt(r2) - 1) : (1 - Math.sqrt(2 - r2))
					d = cx * (((sx + 0.5 + dx) / 2 + x) / w - 0.5) + \
						cy * (((sy + 0.5 + dy) / 2 + y) / h - 0.5) + cam.d

					# Camera rays are pushed forward to start in interior
					camRay = Ray.new(cam.o + d * 140, d.norm)

					# Accumuate radiance
					r = r + radiance(camRay, 0) * (1.0 / samps)
				end

				# Convert radiance to color
				c[i] = c[i] + Vec.new(clamp(r.x), clamp(r.y), clamp(r.z)) * 0.25
			end
		end
	end
end
				

print "\n#{Time.now - start} sec"

f = File.new("image.ppm", 'w')
f.print "P3\n#{w} #{h}\n255\n"
for i in 0..w * h - 1
	f.print "#{toInt(c[i].x)} #{toInt(c[i].y)} #{toInt(c[i].z)}\n"
end
