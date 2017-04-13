function RandomLCG(seed)
	return function ()
		seed = (214013 * seed + 2531011) % 4294967296
		return seed * (1.0 / 4294967296.0)
	end
end

---------------------------------------
Vec = {}
Vec.__index = Vec

function Vec.new(x_, y_, z_)
	local self = { x = x_, y = y_, z = z_}
	setmetatable(self, Vec)
	return self
end
	
function Vec.__add(a, b)
	return Vec.new(a.x + b.x, a.y + b.y, a.z + b.z)
end

function Vec.__sub(a, b)
	return Vec.new(a.x - b.x, a.y - b.y, a.z - b.z)
end

function Vec.__mul(a, b)
	return Vec.new(a.x * b, a.y * b, a.z * b)
end

-- component-wise multiplication
function Vec:mult(b)
	return Vec.new(self.x * b.x, self.y * b.y, self.z * b.z)
end

function Vec:norm()
    return self * (1.0 / math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z))
end

function Vec:dot(b)
	return self.x * b.x + self.y * b.y + self.z * b.z
end

-- cross product
function Vec.__mod(a, b)
    return Vec.new(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
end

Vec.Zero = Vec.new(0, 0, 0)
Vec.XAxis = Vec.new(1, 0, 0)
Vec.YAxis = Vec.new(0, 1, 0)
Vec.ZAxis = Vec.new(0, 0, 1)

---------------------------------------
Refl =
{
    DIFF = 0,
    SPEC = 1,
    REFR = 2
}

---------------------------------------
Ray = {}
Ray.__index = Ray

function Ray.new(o_, d_)
	local self = { o = o_, d = d_ }
	setmetatable(self, Ray)
	return self
end

---------------------------------------
Sphere = {}
Sphere.__index = Sphere

function Sphere.new(rad_, p_, e_, c_, refl_)
	local self = { rad = rad_, p = p_, e = e_, c = c_, refl = refl_ }
	self.sqRad = rad_ * rad_
	self.maxC = math.max(math.max(c_.x, c_.y), c_.z)
	self.cc = c_ * (1.0 / self.maxC)

	setmetatable(self, Sphere)
	return self
end

function Sphere:intersect(r)
	-- Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
	local op = self.p - r.o
	local b = op:dot(r.d)
	local det = b * b - op:dot(op) + self.sqRad
	local eps = 1e-4

	if det < 0 then
		return 0
	else
		local dets = math.sqrt(det)

		if b - dets > eps then
			return b - dets
		elseif b + dets > eps then
			return b + dets
		else
			return 0
		end
	end
end

---------------------------------------
-- Scene: radius, position, emission, color, material
spheres =
{
	Sphere.new(1e5,  Vec.new( 1e5+1,40.8,81.6),  Vec.Zero, Vec.new(.75,.25,.25), Refl.DIFF), --Left
	Sphere.new(1e5,  Vec.new(-1e5+99,40.8,81.6), Vec.Zero, Vec.new(.25,.25,.75), Refl.DIFF), --Rght
	Sphere.new(1e5,  Vec.new(50,40.8, 1e5),      Vec.Zero, Vec.new(.75,.75,.75), Refl.DIFF), --Back
	Sphere.new(1e5,  Vec.new(50,40.8,-1e5+170),  Vec.Zero, Vec.Zero,             Refl.DIFF), --Frnt
	Sphere.new(1e5,  Vec.new(50, 1e5, 81.6),     Vec.Zero, Vec.new(.75,.75,.75), Refl.DIFF), --Botm
	Sphere.new(1e5,  Vec.new(50,-1e5+81.6,81.6), Vec.Zero, Vec.new(.75,.75,.75), Refl.DIFF), --Top
	Sphere.new(16.5, Vec.new(27,16.5,47),        Vec.Zero, Vec.new(1,1,1)*.999,  Refl.SPEC), --Mirr
	Sphere.new(16.5, Vec.new(73,16.5,78),        Vec.Zero, Vec.new(1,1,1)*.999,  Refl.REFR), --Glas
	Sphere.new(600,  Vec.new(50,681.6-.27,81.6), Vec.new(12,12,12), Vec.Zero,    Refl.DIFF)  --Lite
}

rand = RandomLCG(0)

function clamp(x)
	if x < 0 then 
		return 0 
	elseif x > 1 then 
		return 1
	else
		return x 
	end
end

function toInt(x)
	return (clamp(x) ^ (1 / 2.2)) * 255 + .5
end

function intersect(r)
	local t = 1e20
	local obj
	for i, s in ipairs(spheres) do
		local d = s:intersect(r)
		if d ~= 0 and d < t then
			t = d
			obj = s
		end
	end
	return obj, t
end

function radiance(r, depth)
	local obj, t
	obj, t = intersect(r)
	
	if obj == nil then
		return Vec.Zero
	else
		local newDepth = depth + 1
		local isMaxDepth = newDepth > 100

		-- Russian roulette for path termination
		local isUseRR = newDepth > 5
		local isRR = isUseRR and rand() < obj.maxC

		if isMaxDepth or (isUseRR and not isRR) then
			return obj.e
		else
			local f = (isUseRR and isRR) and obj.cc or obj.c
			local x = r.o + r.d * t
			local n = (x - obj.p):norm()
			local nl = (n:dot(r.d) < 0) and n or (n * -1)

            if obj.refl == Refl.DIFF then			-- Ideal DIFFUSE reflection
	            
                local r1 = 2 * math.pi * rand()
                local r2 = rand()
                local r2s = math.sqrt(r2)

                local w = nl
                local wo = (math.abs(w.x) > .1) and Vec.YAxis or Vec.XAxis
                local u = (wo % w):norm()
                local v = w % u

                local d = (u * math.cos(r1) * r2s + v * math.sin(r1) * r2s + w * math.sqrt(1 - r2)):norm()

                return obj.e + f:mult(radiance(Ray.new(x, d), newDepth))                
                
            elseif obj.refl == Refl.SPEC then		-- Ideal SPECULAR reflection
                return obj.e + f:mult(radiance(Ray.new(x, r.d - n * 2 * n:dot(r.d)), newDepth))
                
            else									-- Ideal dielectric REFRACTION
                local reflRay = Ray.new(x, r.d - n * (2 * n:dot(r.d)))
                local into = n:dot(nl) > 0  -- Ray from outside going in?
                local nc = 1
                local nt = 1.5
                local nnt = into and (nc / nt) or (nt / nc)
                local ddn = r.d:dot(nl)
                local cos2t = 1 - nnt * nnt * (1 - ddn * ddn)

                if cos2t < 0 then  -- Total internal reflection
                    return obj.e + f:mult(radiance(reflRay, newDepth))
                else
                    local tdir = (r.d * nnt - n * ((into and 1 or -1) * (ddn * nnt + math.sqrt(cos2t)))):norm()
                    local a = nt - nc
                    local b = nt + nc
                    local R0 = (a * a) / (b * b)
                    local c = 1 - (into and -ddn or tdir:dot(n))
                    local Re = R0 + (1 - R0) * c * c * c * c * c
                    local Tr = 1 - Re
                    local P = .25 + .5 * Re
                    local RP = Re / P
                    local TP = Tr / (1 - P)

                    local result
                    if newDepth > 2 then
                        -- Russian roulette and splitting for selecting reflection and/or refraction
                        if rand() < P then
                            result = radiance(reflRay, newDepth) * RP
                        else
                            result = radiance(Ray.new(x, tdir), newDepth) * TP
                        end
                    else
                        result = radiance(reflRay, newDepth) * Re + radiance(Ray.new(x, tdir), newDepth) * Tr
					end
					
                    return obj.e + f:mult(result)
                end
			end
		end
	end
end

local start = os.clock()

local w = 256
local h = 256

local samps = 25

-- cam pos, dir
local cam = Ray.new(Vec.new(50, 52, 295.6), Vec.new(0, -0.042612, -1):norm())
local cx = Vec.new(w * .5135 / h, 0, 0)
local cy = (cx % cam.d):norm() * .5135

-- final color buffer
local c = {}

-- Loop over image rows
for y = 0, h - 1 do
    io.stderr:write(string.format("\rRendering (%d spp) %5.2f%%", samps * 4, 100 * y / (h - 1)))

    -- Loop cols
    for x = 0, w - 1 do
        local i = (h - y - 1) * w + x
        c[i] = Vec.Zero

        -- 2x2 subpixel rows
        for sy = 0, 1 do
            -- 2x2 subpixel cols
            for sx = 0, 1 do
                local r = Vec.Zero
                for s = 1, samps do
                    local r1 = 2 * rand()
                    local r2 = 2 * rand()
                    local dx = (r1 < 1) and (math.sqrt(r1) - 1) or (1 - math.sqrt(2 - r1))
                    local dy = (r2 < 1) and (math.sqrt(r2) - 1) or (1 - math.sqrt(2 - r2))
                    local d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                              cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d

                    -- Camera rays are pushed forward to start in interior
                    local camRay = Ray.new(cam.o + d * 140, d:norm())

                    -- Accumuate radiance
                    r = r + radiance(camRay, 0) * (1.0 / samps)
                end

                -- Convert radiance to color
                c[i] = c[i] + Vec.new(clamp(r.x), clamp(r.y), clamp(r.z)) * .25
            end
        end
    end
end

print(string.format("\n%f sec", os.clock() - start))

local f = io.open("image.ppm", "w")
f:write(string.format("P3\n%d %d\n%d\n", w, h, 255))
for i = 0, w * h -1 do
	f:write(string.format("%d %d %d\n", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z)))
end
