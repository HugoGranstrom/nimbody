# Port of https://github.com/aprell/tasking-2.0/blob/master/test/nbody3.c
#[
/*
 * The Great Computer Language Shootout
 * http://shootout.alioth.debian.org/
 *
 * contributed by Christoph Bauer
 *
 * adapted by Andreas Prell
 *
 */
]#
import std/[math, times, monotimes, strformat]

type
    Planet = object
        x, y, z: float
        vx, vy, vz: float
        mass: float

const solar_mass = 4*PI*PI
const days_per_year = 365.24
const softening = 0.01


proc advance(bodies, bodies2: ptr UncheckedArray[Planet], i: int, nbodies: int, tick: int, dt: float) =
    var froml, to: ptr UncheckedArray[Planet]
    if tick == 1:
        froml = bodies2
        to = bodies
    else:
        froml = bodies
        to = bodies2
    # copy the internal object to b? Does `=` do that?
    var b = froml[i] #memcpy?
    for j in 0 ..< nbodies:
        let b2 = froml[j]
        if j == i:
            continue
        let dx = b.x - b2.x
        let dy = b.y - b2.y
        let dz = b.z - b2.z
        let dist = sqrt(dx * dx + dy * dy + dz * dz + softening)
        let mag = dt / (dist * dist * dist)
        b.vx -= dx * b2.mass * mag
        b.vy -= dy * b2.mass * mag
        b.vz -= dz * b2.mass * mag
    b.x += dt * b.vx
    b.y += dt * b.vy
    b.z += dt * b.vz
    to[i] = b

proc energy(nbodies: int, bodies: seq[Planet]): float =
    var e: float = 0.0
    for i in 0 ..< nbodies:
        let b = bodies[i]
        e += 0.5 * b.mass * (b.vx*b.vx + b.vy*b.vy + b.vz*b.vz)
        for j in i+1 ..< nbodies:
            let b2 = bodies[j]
            let dx = b.x - b2.x
            let dy = b.y - b2.y
            let dz = b.z - b2.z
            let dist = sqrt(dx*dx + dy*dy + dz*dz)
            e -= (b.mass * b2.mass) / dist
    result = e

proc offset_momentum(nbodies: int, bodies: var seq[Planet]) =
    var px, py, pz: float
    for i in 0 ..< nbodies:
        px += bodies[i].vx * bodies[i].mass
        py += bodies[i].vy * bodies[i].mass
        pz += bodies[i].vz * bodies[i].mass
    bodies[0].vx = - px / solar_mass;
    bodies[0].vy = - py / solar_mass;
    bodies[0].vz = - pz / solar_mass;

proc fscanf(f: File, s: cstring)
  {.varargs,
    importc: "fscanf",
    header: "<stdio.h>".}

proc simul_mono(n, N: int, bodies, bodies2: var seq[Planet], dt: float64) =

  let bodies_ptr = cast[ptr UncheckedArray[Planet]](bodies[0].unsafeAddr)
  let bodies2_ptr = cast[ptr UncheckedArray[Planet]](bodies2[0].unsafeAddr)

  echo "Simulating ", N, " bodies for ", n, " timesteps"
  offset_momentum(N, bodies)
  echo energy(N, bodies)
  let startTime = getMonotime()
  for i in 0 ..< n:
      for j in 0 ..< N:
          advance(bodies_ptr, bodies2_ptr, j, N, i mod 2, dt)
  let endTime = getMonotime()
  let elapsed = inMilliseconds(endTime - startTime)
  echo &"Elapsed wall time (single-threaded): {elapsed:>6} ms ({elapsed.float64/n.float64:>6.3f} ms per iteration)"

when compileOption("threads"):
  import weave, cpuinfo

  proc simul_weave(n, N: int, bodies, bodies2: var seq[Planet], dt: float64) =

    let bodies_ptr = cast[ptr UncheckedArray[Planet]](bodies[0].unsafeAddr)
    let bodies2_ptr = cast[ptr UncheckedArray[Planet]](bodies2[0].unsafeAddr)

    init(Weave)
    echo "Simulating ", N, " bodies for ", n, " timesteps"
    offset_momentum(N, bodies)
    echo energy(N, bodies)
    let startTime = getMonotime()
    for i in 0 ..< n:
        parallelFor j in 0 ..< N:
            captures: {bodies_ptr, bodies2_ptr, N, i, dt}
            advance(bodies_ptr, bodies2_ptr, j, N, i mod 2, dt)
        syncRoot(Weave)
    let endTime = getMonotime()
    let elapsed = inMilliseconds(endTime - startTime)
    echo &"Elapsed wall time (Weave): {elapsed:>6} ms ({elapsed.float64/n.float64:>6.3f} ms per iteration)"

proc main() =
    let dt = 0.001
    let n = 1000 # number of iterations
    var input: File
    input = open("bodies.txt")
    var N_c: cint
    fscanf(input, "%d", addr N_c)
    let N = int(N_c)
    var bodies = newSeq[Planet](N)
    for i in 0 ..< N:
        var a: cdouble
        fscanf(input, "%lf", addr a)
        bodies[i].mass = float(a)
    for i in 0 ..< N:
        var a: cdouble
        fscanf(input, "%lf", addr a)
        bodies[i].x = float(a)
    for i in 0 ..< N:
        var a: cdouble
        fscanf(input, "%lf", addr a)
        bodies[i].y = float(a)
    for i in 0 ..< N:
        var a: cdouble
        fscanf(input, "%lf", addr a)
        bodies[i].z = float(a)
    for i in 0 ..< N:
        var a: cdouble
        fscanf(input, "%lf", addr a)
        bodies[i].vx = float(a)
    for i in 0 ..< N:
        var a: cdouble
        fscanf(input, "%lf", addr a)
        bodies[i].vy = float(a)
    for i in 0 ..< N:
        var a: cdouble
        fscanf(input, "%lf", addr a)
        bodies[i].vz = float(a)
    input.close()
    var bodies2 = bodies # copy?

    simul_mono(n, N, bodies, bodies2, dt)
    simul_weave(n, N, bodies, bodies2, dt)

    echo energy(N, bodies)

main()
