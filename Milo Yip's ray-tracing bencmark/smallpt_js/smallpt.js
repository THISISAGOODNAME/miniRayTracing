﻿function RandomLCG(seed) {
    return function() {
        seed = (214013 * seed + 2531011) % 0x100000000;
        return seed * (1.0 / 4294967296.0);
    };
}

function Vec(x, y, z) { this.x = x; this.y = y; this.z = z; }

Vec.prototype =
{
    add: function(b) {
        return new Vec(this.x + b.x, this.y + b.y, this.z + b.z);
    },

    sub: function(b) {
        return new Vec(this.x - b.x, this.y - b.y, this.z - b.z);
    },

    mul: function(b) {
        return new Vec(this.x * b, this.y * b, this.z * b);
    },

    mult: function(b) {
        return new Vec(this.x * b.x, this.y * b.y, this.z * b.z);
    },

    norm: function() {
        return this.mul(1.0 / Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z));
    },

    dot: function(b) {
        return this.x * b.x + this.y * b.y + this.z * b.z;
    },

    cross: function(b) {
        return new Vec(this.y * b.z - this.z * b.y, this.z * b.x - this.x * b.z, this.x * b.y - this.y * b.x);
    }
};

Vec.Zero = new Vec(0, 0, 0)
Vec.XAxis = new Vec(1, 0, 0)
Vec.YAxis = new Vec(0, 1, 0)
Vec.ZAxis = new Vec(0, 0, 1)

Refl = {
    DIFF: 0,
    SPEC: 1,
    REFR: 2
};

function Ray(o, d) { this.o = o; this.d = d }

function Sphere(rad, p, e, c, refl) {
    this.rad = rad;
    this.p = p;
    this.e = e;
    this.c = c;
    this.refl = refl;

    this.sqRad = rad * rad;
    this.maxC = Math.max(Math.max(c.x, c.y), c.z);
    this.cc = c.mul(1.0 / this.maxC);
}

Sphere.prototype =
{
    intersect: function(r) {
        // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        var op = this.p.sub(r.o);
        var b = op.dot(r.d);
        var det = b * b - op.dot(op) + this.sqRad;
        var eps = 1e-4;

        if (det < 0)
            return 0;
        else {
            var dets = Math.sqrt(det)

            if (b - dets > eps)
                return b - dets;
            else if (b + dets > eps)
                return b + dets;
            else
                return 0;
        }
    }
};

// Scene: radius, position, emission, color, material
spheres = [
    new Sphere(1e5, new Vec(1e5 + 1, 40.8, 81.6), Vec.Zero, new Vec(.75, .25, .25), Refl.DIFF), //Left
    new Sphere(1e5, new Vec(-1e5 + 99, 40.8, 81.6), Vec.Zero, new Vec(.25, .25, .75), Refl.DIFF), //Rght
    new Sphere(1e5, new Vec(50, 40.8, 1e5), Vec.Zero, new Vec(.75, .75, .75), Refl.DIFF), //Back
    new Sphere(1e5, new Vec(50, 40.8, -1e5 + 170), Vec.Zero, Vec.Zero, Refl.DIFF), //Frnt
    new Sphere(1e5, new Vec(50, 1e5, 81.6), Vec.Zero, new Vec(.75, .75, .75), Refl.DIFF), //Botm
    new Sphere(1e5, new Vec(50, -1e5 + 81.6, 81.6), Vec.Zero, new Vec(.75, .75, .75), Refl.DIFF), //Top
    new Sphere(16.5, new Vec(27, 16.5, 47), Vec.Zero, new Vec(1, 1, 1).mul(.999), Refl.SPEC), //Mirr
    new Sphere(16.5, new Vec(73, 16.5, 78), Vec.Zero, new Vec(1, 1, 1).mul(.999), Refl.REFR), //Glas
    new Sphere(600, new Vec(50, 681.6 - .27, 81.6), new Vec(12, 12, 12), Vec.Zero, Refl.DIFF)  //Lite
];

var rand = RandomLCG(0)

function clamp(x) {
    if (x < 0)
        return 0;
    else if (x > 1)
        return 1;
    else
        return x;
}

function toInt(x) {
    return Math.pow(clamp(x), 1 / 2.2) * 255 + .5
}

function intersect(r) {
    var t = 1e20;
    var obj;

    for (var i in spheres) {
        var s = spheres[i];
        var d = s.intersect(r);
        if (d != 0 && d < t) {
            t = d;
            obj = s;
        }
    }
    return [obj, t];
}

function radiance(r, depth) {
    var ires = intersect(r);
    var obj = ires[0];
    var t = ires[1];   // distance to intersection

    if (obj == null) {
        return Vec.Zero;       // if miss, return black
    }
    else {
        var newDepth = depth + 1;
        var isMaxDepth = newDepth > 100;

        // Russian roulette for path termination
        var isUseRR = newDepth > 5;
        var isRR = isUseRR && rand() < obj.maxC;

        if (isMaxDepth || (isUseRR && !isRR))
            return obj.e;
        else {
            var f = (isUseRR && isRR) ? obj.cc : obj.c;
            var x = r.o.add(r.d.mul(t));
            var n = x.sub(obj.p).norm();
            var nl = n.dot(r.d) < 0 ? n : n.mul(-1);

            if (obj.refl == Refl.DIFF) // Ideal DIFFUSE reflection
            {
                var r1 = 2 * Math.PI * rand();
                var r2 = rand();
                var r2s = Math.sqrt(r2);

                var w = nl;
                var wo = Math.abs(w.x) > .1 ? new Vec(0, 1, 0) : new Vec(1, 0, 0);
                var u = wo.cross(w).norm();
                var v = w.cross(u);

                var d = (u.mul(Math.cos(r1) * r2s).add(v.mul(Math.sin(r1) * r2s)).add(w.mul(Math.sqrt(1 - r2)))).norm();

                return obj.e.add(f.mult(radiance(new Ray(x, d), newDepth)));
            }
            else if (obj.refl == Refl.SPEC) // Ideal SPECULAR reflection
            {
                return obj.e.add(f.mult(radiance(new Ray(x, r.d.sub(n.mul(2 * n.dot(r.d)))), newDepth)));
            }
            else // Ideal dielectric REFRACTION
            {
                var reflRay = new Ray(x, r.d.sub(n.mul(2 * n.dot(r.d))));
                var into = n.dot(nl) > 0; // var from outside going in?
                var nc = 1;
                var nt = 1.5;
                var nnt = into ? nc / nt : nt / nc;
                var ddn = r.d.dot(nl);
                var cos2t = 1 - nnt * nnt * (1 - ddn * ddn);

                if (cos2t < 0)	// Total internal reflection
                {
                    return obj.e.add(f.mult(radiance(reflRay, newDepth)));
                }
                else {
                    var tdir = (r.d.mul(nnt).sub(n.mul((into ? 1 : -1) * (ddn * nnt + Math.sqrt(cos2t))))).norm();
                    var a = nt - nc;
                    var b = nt + nc;
                    var R0 = a * a / (b * b);
                    var c = 1 - (into ? -ddn : tdir.dot(n));
                    var Re = R0 + (1 - R0) * c * c * c * c * c;
                    var Tr = 1 - Re;
                    var P = .25 + .5 * Re;
                    var RP = Re / P;
                    var TP = Tr / (1 - P);

                    var result;
                    if (newDepth > 2) {
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

function render(canvas, status) {
    var start = new Date().getTime();
    var w = canvas.attributes.width.value;
    var h = canvas.attributes.height.value;
    var samps = 25;

    // cam pos, dir
    var cam = new Ray(new Vec(50, 52, 295.6), new Vec(0, -0.042612, -1).norm());
    var cx = new Vec(w * .5135 / h, 0, 0);
    var cy = (cx.cross(cam.d)).norm().mul(.5135);

    // final color buffer
    var c = new Array(w * h);
    for (var i = 0; i < w * h; i++)
        c[i] = Vec.Zero;

    // Output
    var ctx = canvas.getContext("2d");
    var imgdata = ctx.getImageData(0, 0, w, h);
    var pixels = imgdata.data;

    // Loop over image rows
    var y = 0;
    setTimeout(renderLine, 0);
    
    function renderLine()
    {
        status.innerHTML = "Rendering (" + samps * 4 + " spp) " + (100.0 * y / (h - 1)).toFixed(2) + "%";
        
        // Loop cols
        for (var x = 0; x < w; x++) {
            // 2x2 subpixel rows
            for (var sy = 0; sy < 2; sy++) {
                var i = (h - y - 1) * w + x;

                // 2x2 subpixel cols
                for (var sx = 0; sx < 2; sx++) {
                    var r = Vec.Zero;
                    for (var s = 0; s < samps; s++) {
                        var r1 = 2 * rand();
                        var r2 = 2 * rand();
                        var dx = r1 < 1 ? Math.sqrt(r1) - 1 : 1 - Math.sqrt(2 - r1);
                        var dy = r2 < 1 ? Math.sqrt(r2) - 1 : 1 - Math.sqrt(2 - r2);

                        var d = cx.mul(((sx + .5 + dx) / 2 + x) / w - .5).add(
					    cy.mul(((sy + .5 + dy) / 2 + y) / h - .5)).add(cam.d);

                        // Camera rays are pushed forward to start in interior
                        var camRay = new Ray(cam.o.add(d.mul(140)), d.norm());

                        // Accumuate radiance
                        r = r.add(radiance(camRay, 0).mul(1.0 / samps));
                    }

                    // Convert radiance to color
                    c[i] = c[i].add((new Vec(clamp(r.x), clamp(r.y), clamp(r.z))).mul(.25));
                }
            }
        }

        renderOutput();

        y++;
        if (y < h)
            setTimeout(renderLine, 0);
        else
            status.innerHTML = (new Date().getTime() - start) / 1000 + " sec";
    }

    function renderOutput() {
        var i = (h - y - 1) * w * 4, j = (h - y - 1) * w;
        for (var x = 0; x < w; x++) {
            pixels[i++] = toInt(c[j].x);
            pixels[i++] = toInt(c[j].y);
            pixels[i++] = toInt(c[j].z);
            pixels[i++] = 255;
            j++;
        }

        ctx.putImageData(imgdata, 0, 0);
    }
}