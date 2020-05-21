import math
import arraymancer, weave, numericalnim

proc normalMain(posStart, velStart: Tensor[float]): Tensor[float] =
    let N = posStart.shape[0]
    let M = posStart.shape[1]
    var pos = posStart.clone()
    var vel = velStart.clone()
    var acc = zeros[float]([N, M])
    var t = 0.0
    let dt: float = 1e-2
    let G = 1e-1
    let T = 1.0
    while t < T:
        acc *= 0.0
        # Calculate accelerations on all bodies
        for i in 0 ..< N:
            for j in 0 ..< N:
                if i != j:
                    let r_vec = pos[j,_] - pos[i,_]
                    let r3 = pow(sum(r_vec *. r_vec), 3/2)
                    for k in 0 ..< 3:
                        acc[i,k] += r_vec[0, k] * (G / r3)
        vel += dt * acc
        pos += dt * vel
        t += dt
    result = pos

proc normalMain2(posStart, velStart: Tensor[float]): Tensor[float] =
    let N = posStart.shape[0]
    let M = posStart.shape[1]
    var pos = posStart.clone()
    var vel = velStart.clone()
    var acc = zeros[float]([N, M])
    let p = pos.dataArray()
    let v = vel.dataArray()
    let a = acc.dataArray()
    var t = 0.0
    let dt: float = 1e-2
    let G = 1e-1
    let T = 1.0
    while t < T:
        acc *= 0.0
        # Calculate acceleration from j:th body on i:th body.
        for i in 0 ..< N:
            for j in 0 ..< N:
                if i != j:
                    let diffx = p[j*M] - p[i*M]
                    let diffy = p[j*M+1] - p[i*M+1]
                    let diffz = p[j*M+2] - p[i*M+2]
                    let r3 = pow(diffx*diffx + diffy*diffy + diffz*diffz, 3/2)
                    let factor = G/r3
                    a[i*M] += factor * diffx
                    a[i*M+1] += factor * diffy
                    a[i*M+2] += factor * diffz
        vel += dt * acc
        pos += dt * vel
        t += dt
    result = pos


proc parallelMain(posStart, velStart: Tensor[float]): Tensor[float] =
    let N = posStart.shape[0]
    let M = posStart.shape[1]
    var pos = posStart.clone()
    var vel = velStart.clone()
    var acc = zeros[float]([N, M])
    let p = pos.dataArray()
    let v = vel.dataArray()
    let a = acc.dataArray()
    var t = 0.0
    let dt: float = 1e-2
    let G = 1e-1
    let T = 1.0
    init(Weave)
    while t < T:
        acc *= 0.0
        # Calculate acceleration from j:th body on i:th body.
        parallelFor i in 0 ..< N:
            captures: {M, N, p, a, G}
            for j in 0 ..< N:
            #parallelFor j in 0 ..< N:
                #captures: {i, M, N, p, a, G}
                if i != j:
                    let diffx = p[j*M] - p[i*M]
                    let diffy = p[j*M+1] - p[i*M+1]
                    let diffz = p[j*M+2] - p[i*M+2]
                    let r3 = pow(diffx*diffx + diffy*diffy + diffz*diffz, 3/2)
                    let factor = G/r3
                    a[i*M] += factor * diffx
                    a[i*M+1] += factor * diffy
                    a[i*M+2] += factor * diffz
        syncRoot(Weave)
        vel += dt * acc
        pos += dt * vel
        t += dt
    exit(Weave)
    result = pos

let N = 1000
let p = randomTensor[float]([N, 3], 10) -. 5.0
let v = randomTensor[float]([N, 3], 2) -. 1
#let posInput = @[@[0.0, 0.0, 0.0], @[1.0, 1.0, 1.0], @[-1.0, 1.0, -1.0]].toTensor()
#let velInput = @[@[0.0, 0.0, 0.0], @[-1.0, 1.0, 0.0], @[0.0, -1.0, 0.0]].toTensor()
#echo parallelMain(p, v)[N-3.._,_]
#echo normalMain2(p, v)[N-3.._,_]
#echo normalMain(p, v)[N-3.._,_]
timeit(parallelMain(p, v), 5, "Parallel")
timeit(normalMain(p, v), 3,   "Normal  ")
timeit(normalMain2(p, v), 5,  "Normal2 ")

