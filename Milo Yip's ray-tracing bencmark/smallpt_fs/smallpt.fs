// Learn more about F# at http://fsharp.net

open System
open System.IO
open System.Diagnostics

let mutable seed = 0u;
let inline randomLCG() =
    seed <- (214013u * seed + 2531011u) % 0xFFFFFFFFu
    (double seed) / 4294967296.0

type Vec = { x: double; y: double; z: double; }
let inline makeVec x y z = { x = x; y = y; z = z; }
let inline (++) (a: Vec) (b: Vec) = makeVec (a.x + b.x) (a.y + b.y) (a.z + b.z)
let inline (--) (a: Vec) (b: Vec) = makeVec (a.x - b.x) (a.y - b.y) (a.z - b.z)
let inline ( ** ) (a: Vec) (b: double) = makeVec (a.x * b) (a.y * b) (a.z * b)
let inline mult (a: Vec) (b: Vec) = makeVec (a.x * b.x) (a.y * b.y) (a.z * b.z)
let inline norm (a: Vec) = a ** (1.0 / sqrt(a.x * a.x + a.y * a.y + a.z * a.z))
let inline dot (a: Vec) (b: Vec) = a.x * b.x + a.y * b.y + a.z * b.z;
let inline (%%) (a: Vec) (b: Vec) = makeVec (a.y * b.z - a.z * b.y) (a.z * b.x - a.x * b.z) (a.x * b.y - a.y * b.x)
let VecZero = makeVec 0.0 0.0 0.0
let VecXAxis = makeVec 1.0 0.0 0.0
let VecYAxis = makeVec 0.0 1.0 0.0
let VecZAxis = makeVec 0.0 0.0 1.0

type Refl_t =
| DIFF = 0x0001 
| SPEC = 0x0002
| REFR = 0x0004

type Ray = { o: Vec; d: Vec; }
let inline makeRay o d = { o = o; d = d }

type Sphere =
    val public rad: double
    val public p: Vec
    val public e: Vec
    val public c: Vec
    val public refl: Refl_t
    val public maxC: double
    val public cc: Vec
    val public sqRad: double

    new (rad_: double, p_: Vec, e_: Vec, c_: Vec, refl_: Refl_t) =
        let maxC_ = max c_.x c_.y |> max c_.z
        { rad = rad_; p = p_; e = e_; c = c_; refl = refl_;
          sqRad = rad_ * rad_; maxC = maxC_; cc = c_ ** (1.0 / maxC_) }

    member inline s.intersect(r: Ray) =
        // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        let op = s.p -- r.o
        let eps = 1e-4
        let b = dot op r.d
        let det = (b * b) - (dot op op) + (s.rad * s.rad)
        
        if det < 0.0 then
            0.0
        else
            let dets = sqrt(det)
            if b - dets > eps then 
                b - dets
            elif b + dets > eps then
                b + dets;
            else
                0.0
let inline makeSphere rad p e c refl = new Sphere(rad, p, e, c, refl)

[<EntryPoint>]
let main args = 
    let inline rand() = randomLCG()

    let spheres = [|
        makeSphere 1e5 (makeVec (1e5 + 1.0) 40.8 81.6) VecZero (makeVec 0.75 0.25 0.25) Refl_t.DIFF; //Left
        makeSphere 1e5 (makeVec (-1e5 + 99.0) 40.8 81.6) VecZero (makeVec 0.25 0.25 0.75) Refl_t.DIFF; //Rght
        makeSphere 1e5 (makeVec 50.0 40.8 1e5) VecZero (makeVec 0.75 0.75 0.75) Refl_t.DIFF; //Back
        makeSphere 1e5 (makeVec 50.0 40.8 (-1e5 + 170.0)) VecZero VecZero Refl_t.DIFF; //Frnt
        makeSphere 1e5 (makeVec 50.0 1e5 81.6) VecZero (makeVec 0.75 0.75 0.75) Refl_t.DIFF; //Botm
        makeSphere 1e5 (makeVec 50.0 (-1e5 + 81.6) 81.6) VecZero (makeVec 0.75 0.75 0.75) Refl_t.DIFF; //Top
        makeSphere 16.5 (makeVec 27.0 16.5 47.0) VecZero ((makeVec 1.0 1.0 1.0) ** 0.999) Refl_t.SPEC; //Mirr
        makeSphere 16.5 (makeVec 73.0 16.5 78.0) VecZero ((makeVec 1.0 1.0 1.0) ** 0.999) Refl_t.REFR; //Glas
        makeSphere 600.0 (makeVec 50.0 (681.6 - 0.27) 81.6) (makeVec 12.0 12.0 12.0) VecZero Refl_t.DIFF; //Lite
    |]

    let inline clamp x =
        if x < 0.0 then 0.0
        elif x > 1.0 then 1.0 else x

    let inline toInt x = int (Math.Pow(clamp x, 1.0 / 2.2) * 255.0 + 0.5)

    let inline intersect (r: Ray) =
        let mutable t = 1e20
        let mutable inf = 1e20
        let mutable ret = None

        for s in spheres do
            let d = s.intersect(r)
            if (d <> 0.0 && d < t) then
                t <- d
                ret <- Some s
        (ret, t)

    let rec radiance (r: Ray) depth = 
        // distance to intersection
        let optionObj, t = intersect r
        match optionObj with
        | None -> VecZero // if miss, return black
        | Some obj ->
            let newDepth = depth + 1
            let isMaxDepth = newDepth > 100
            
            // Russian roulette for path termination
            let isUseRR = newDepth > 5
            let isRR = isUseRR && rand() < obj.maxC

            if isMaxDepth || (isUseRR && (not isRR)) then
                obj.e
            else
                let f = if (isUseRR && isRR) then obj.cc else obj.c
                let x = r.o ++ (r.d ** t);
                let n = norm (x -- obj.p)
                let nl = if dot n r.d < 0.0 then n else n ** -1.0

                if obj.refl = Refl_t.DIFF then // Ideal DIFFUSE reflection
                    let r1 = 2.0 * Math.PI * rand()
                    let r2 = rand()
                    let r2s = sqrt r2

                    let w = nl;
                    let wo = if abs w.x > 0.1 then VecYAxis else VecXAxis
                    let u = norm (wo %% w)
                    let v = w %% u

                    let d = norm ((u ** (cos r1)) ** r2s ++ (v ** (sin r1)) ** r2s ++ w ** (sqrt (1.0 - r2)))
                    obj.e ++ mult f (radiance (makeRay x d) newDepth)
                elif obj.refl = Refl_t.SPEC then // Ideal SPECULAR reflection
                    obj.e ++ mult f (radiance (makeRay x (r.d -- (n ** 2.0) ** (dot n r.d))) newDepth)
                else // Ideal dielectric REFRACTION
                    let reflRay = makeRay x (r.d -- (n ** 2.0) ** (dot n r.d))
                    let into = dot n nl > 0.0 // Ray from outside going in?
                    let nc = 1.0
                    let nt = 1.5
                    let nnt = if into then nc / nt else nt / nc
                    let ddn = dot r.d nl
                    let cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn)

                    if cos2t < 0.0 then // Total internal reflection
                        obj.e ++ (mult f (radiance reflRay newDepth))
                    else
                        let t = if into then 1.0 else -1.0
                        let tdir = norm (r.d ** nnt -- n ** (t * (ddn * nnt + (sqrt cos2t))))
                        let a = nt - nc
                        let b = nt + nc
                        let R0 = (a * a) / (b * b)
                        let c = 1.0 - if into then -ddn else (dot tdir n)
                        let Re = R0 + (1.0 - R0) * c * c * c * c * c
                        let Tr = 1.0 - Re
                        let P = 0.25 + 0.5 * Re
                        let RP = Re / P
                        let TP = Tr / (1.0 - P)

                        let result = 
                            if newDepth > 2 then
                                // Russian roulette and splitting for selecting reflection and/or refraction
                                if rand() < P then
                                    (radiance reflRay depth) ** RP
                                else
                                    (radiance (makeRay x tdir) depth) ** TP
                            else
                                ((radiance reflRay newDepth) ** Re) ++ ((radiance (makeRay x tdir) newDepth) ** Tr)

                        obj.e ++ (mult f result)

    // main
    let watch = Stopwatch.StartNew()

    let w = 256;
    let h = 256;

    let samps = if args.Length = 2 then Int32.Parse(args.[1]) / 4 else 25; // # samples

    // cam pos
    let cam = makeRay (makeVec 50.0 52.0 295.6) (makeVec 0.0 -0.042612 -1.0 |> norm)

    // cam dir
    let cx = makeVec (double w * 0.5135 / double h) 0.0 0.0
    let cy = ((cx %% cam.d) |> norm) ** 0.5135

    // final color buffer
    let c = Array.init (w * h) (fun _ -> VecZero)
    
    // Loop over image rows
    for y = 0 to h - 1 do
        Console.Write("\rRendering ({0} spp) {1:F2}%", samps * 4 |> box, 100.0 * double y / double (h - 1) |> box)
        // Loop cols
        for x = 0 to w - 1 do
            let i = (h - y - 1) * w + x;
            // 2x2 subpixel rows
            for sy = 0 to 1 do
                // 2x2 subpixel rows
                for sx = 0 to 1 do
                    let mutable r = VecZero
                    for s = 0 to samps - 1 do
                        let r1 = 2.0 * rand()
                        let dx = if r1 < 1.0 then sqrt(r1) - 1.0 else 1.0 - sqrt(2.0 - r1)
                        let r2 = 2.0 * rand()
                        let dy = if r2 < 1.0 then sqrt(r2) - 1.0 else 1.0 - sqrt(2.0 - r2)
                        let d = 
                            cx ** (((double sx + 0.5 + double dx) / 2.0 + double x) / double w - 0.5) ++
                            cy ** (((double sy + 0.5 + double dy) / 2.0 + double y) / double h - 0.5) ++ cam.d;

                        // Camera rays are pushed forward to start in interior
                        let camRay = makeRay (cam.o ++ d ** 140.0) (norm d)

                        // Accumuate radiance
                        r <- r ++ (radiance camRay 0) ** (1.0 / double samps)

                    // Convert radiance to color
                    c.[i] <- c.[i] ++ (makeVec (clamp r.x) (clamp r.y) (clamp r.z)) ** 0.25;

    watch.Elapsed.ToString() |> printfn "\n%s sec"

    using (new StreamWriter("image.ppm")) (fun sw ->
        sw.Write("P3\r\n{0} {1}\r\n{2}\r\n", w, h, 255);
        for i = 0 to w * h - 1 do
            sw.Write("{0} {1} {2}\r\n", toInt(c.[i].x), toInt(c.[i].y), toInt(c.[i].z))
    )

    0