import math, tables, random, std/monotimes, times
import weave
import types, accMethods, integrators

export types, accMethods, integrators



############################
#      Integrators         #
############################

proc kineticEnergy*(bodies: ptr UncheckedArray[Body], nBodies: int64): float64 {.inline.} =
    for i in 0 ..< nBodies:
        let body = bodies[i]
        result += 0.5 * body.mass * dot(body.vel, body.vel)


proc gravityPotentialEnergy*(bodies: ptr UncheckedArray[Body], nBodies: int64, G: float64): float64 {.inline.} =
    result = 0.0
    for i in 0 ..< nBodies:
        let body1 = bodies[i]
        for j in i+1 ..< nBodies:
            let body2 = bodies[j]
            let r = mag(body1.pos - body2.pos)
            result -= G * body1.mass * body2.mass / r    

proc kineticEnergyParallel*(bodies: ptr UncheckedArray[Body], nBodies: int64): float64 {.inline.} =
    var waitableSum: Flowvar[float64]
    parallelFor i in 0 .. nBodies.int:
        captures: {bodies}
        reduce(waitableSum):
            prologue:
                var localSum: float64 = 0
            fold:
                let body = bodies[i]
                localSum += 0.5 * body.mass * dot(body.vel, body.vel)
            merge(remoteSum):
                localSum += sync(remoteSum)
            return localSum 
    result = sync(waitableSum)


proc gravityPotentialEnergyParallel*(bodies: ptr UncheckedArray[Body], nBodies: int64, G: float64): float64 {.inline.} =
    var waitableSum: Flowvar[float64]
    parallelFor i in 0 .. nBodies.int:
        captures: {bodies, nBodies, G}
        reduce(waitableSum):
            prologue:
                var localSum: float64 = 0
            fold:
                let body1 = bodies[i]
                for j in i+1 ..< nBodies:
                    let body2 = bodies[j]
                    let r = mag(body1.pos - body2.pos)
                    localSum -= G * body1.mass * body2.mass / r
            merge(remoteSum):
                localSum += sync(remoteSum)
            return localSum 
    result = sync(waitableSum)


proc nbodySimulation*[T, U, V](initBodies: openArray[Body], calcAcc: CalcAccProc[T, U, V], integrator: IntegratorProc[T, U, V], simul_time, dt: float64, customMem1: PersistantMemory[T] = nil, customMem2: PersistantMemory[U] = nil, customMem3: PersistantMemory[V] = nil): (seq[Body], PersistantMemory[float], PersistantMemory[T], PersistantMemory[U], PersistantMemory[V]) =
    # what other parameters do we need to pass?
    let floatMem = PersistantMemory[float64]()
    var customMem1Internal = customMem1
    var customMem2Internal = customMem2
    var customMem3Internal = customMem3
    if isNil(customMem1Internal):
        customMem1Internal = PersistantMemory[T]()
    if isNil(customMem2Internal):
        customMem2Internal = PersistantMemory[U]()
    if isNil(customMem3Internal):
        customMem3Internal = PersistantMemory[V]()
    var bodies = @initBodies # make a copy and convert to seq
    let nBodies = bodies.len
    let bodies_ptr = cast[ptr UncheckedArray[Body]](bodies[0].unsafeaddr) # make an uncheckedarray out of it
    let n_whole_steps = abs(floor(simul_time / dt)).int64 # do n full steps and then a correction-step
    #let startTime = getMonoTime()
    #var kin = kineticEnergyParallel(bodies_ptr, nBodies)
    #var pot = gravityPotentialEnergyParallel(bodies_ptr, nBodies, 1e-3)
    #var tot = kin + pot
    #echo "Befor: Kinetic = ", kin, " Potential = ", pot, " Total = ", kin+pot
    for timestep in 0 ..< n_whole_steps:
        integrator(timestep, bodies_ptr, nBodies, calcAcc, dt, floatMem, customMem1Internal, customMem2Internal, customMem3Internal)
        # pre force hook
        # integrator_step
        # after force hook (ex: calculate energy, save to hdf5)
    # do correction step to end on simul_time
    let last_dt = simul_time - n_whole_steps.float64 * dt
    integrator(-1, bodies_ptr, bodies.len, calcAcc, last_dt, floatMem, customMem1Internal, customMem2Internal, customMem3Internal)
    # timestep = -1 on the last iteration
    #kin = kineticEnergyParallel(bodies_ptr, nBodies)
    #pot = gravityPotentialEnergyParallel(bodies_ptr, nBodies, 1e-3)
    #echo "After: Kinetic = ", kin, " Potential = ", pot, " Total = ", kin+pot
    #echo "diff: ", tot - (kin + pot)
    #let endTime = getMonoTime()
    #let elapsed = inMilliSeconds(endTime - startTime)
    #echo "Elapsed time: ", elapsed, " ms"
        
    #echo "Energy diff: ", tot - (kin+pot)
    return (bodies, floatMem, customMem1Internal, customMem2Internal, customMem3Internal)
    
when isMainModule:
    var a = Body(mass: 1.0)
    let b = PersistantMemory[float64]()
    let N = 1_000_00
    var bodies = newSeq[Body](N)
    let bodies_ptr = cast[ptr UncheckedArray[Body]](bodies[0].unsafeAddr)
    var totalmass: float64
    var com = Vector3()
    for i in 0 ..< N:
        let x = (rand(1.0) - 0.5) * 20
        let y = (rand(1.0) - 0.5) * 20
        let z = (rand(1.0) - 0.5) * 20
        let pos = Vector3(x: x, y: y, z: z)
        let vel = Vector3(x: -y, y: x, z: z) / 1000
        let mass = rand(1e-2)
        totalmass += mass
        com += pos * mass
        bodies[i] = Body(pos: pos, vel: vel, mass: mass)
    #echo "com: ", com / totalmass
    #echo "Total mass: ", totalmass
    init(Weave)
    let (bodiesFinal1, floatM1, M11, M21, M31) = nbodySimulation[float64, float64, int64](initBodies = bodies, calcAcc = barnesHutAccParallel, integrator = velocityVerletStep, simul_time = 1e0, dt = 0.3, customMem1 = b)
    let (bodiesFinal2, floatM2, M12, M22, M32) = nbodySimulation[float64, float64, int64](initBodies = bodies, calcAcc = barnesHutAccParallel, integrator = velocityVerletStep, simul_time = 1e0, dt = 0.3, customMem1 = b)
    discard nbodySimulation[float64, float64, int64](initBodies = bodies, calcAcc = barnesHutAcc, integrator = velocityVerletStep, simul_time = 1e0, dt = 0.3, customMem1 = b)
    #let (bodiesFinal3, floatM3, M13, M23, M33) = nbodySimulation[float64, float64, int64](initBodies = bodies, calcAcc = barnesHutAcc, integrator = velocityVerletStepParallel, simul_time = 1e0, dt = 1e0, customMem1 = b)
    let b1 = bodiesFinal1[0]
    let b2 = bodiesFinal1[0]
    exit(Weave)
    #echo bodies[0]
    #echo b1
    #echo b2
    echo "b1 == b2: ", b1 == b2
    #echo "b1 == b0: ", b1 == bodies[0]
    
    #echo M2[]
    #echo M2
    #echo bodies
    #echo sizeof(c[0])