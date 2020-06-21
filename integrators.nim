import weave
import types

# Verlet

proc velocityVerletStep*[T, U, V](timestep: int64, bodies: ptr UncheckedArray[Body], nBodies: int64, calcAcc: CalcAccProc[T, U, V], dt: float64, floatMem: PersistantMemory[float64], customMem1: PersistantMemory[T], customMem2: PersistantMemory[U], customMem3: PersistantMemory[V]) =
    if timestep == 0:
        let accFirst = calcAcc(timestep, bodies, nBodies, floatMem, customMem1, customMem2, customMem3)
        for i in 0 ..< nBodies:
            bodies[i].memory.vector = accFirst[i]
    for i in 0 ..< nBodies:
        bodies[i].pos += bodies[i].vel * dt + (0.5 * dt * dt) * bodies[i].memory.vector 
    let acc = calcAcc(timestep, bodies, nBodies, floatMem, customMem1, customMem2, customMem3)
    for i in 0 ..< nBodies:
        bodies[i].vel += (0.5 * dt) * (bodies[i].memory.vector + acc[i])
        bodies[i].memory.vector = acc[i]

proc velocityVerletStepParallel*[T, U, V](timestep: int64, bodies: ptr UncheckedArray[Body], nBodies: int64, calcAcc: CalcAccProc[T, U, V], dt: float64, floatMem: PersistantMemory[float64], customMem1: PersistantMemory[T], customMem2: PersistantMemory[U], customMem3: PersistantMemory[V]) =
    let dt = dt
    if timestep == 0:
        let accFirst = calcAcc(timestep, bodies, nBodies, floatMem, customMem1, customMem2, customMem3)
        let acc_ptr = cast[ptr UncheckedArray[Vector3]](accFirst[0].unsafeAddr)
        syncScope():
            parallelFor i in 0 ..< nBodies.int:
                captures: {bodies, acc_ptr}
                bodies[i].memory.vector = acc_ptr[i]
    syncScope():
        parallelFor i in 0 ..< nBodies.int:
            captures: {bodies, dt}
            bodies[i].pos += bodies[i].vel * dt + 0.5 * bodies[i].memory.vector * dt * dt
    let acc = calcAcc(timestep, bodies, nBodies, floatMem, customMem1, customMem2, customMem3)
    let acc_ptr = cast[ptr UncheckedArray[Vector3]](acc[0].unsafeAddr)
    syncScope():
        parallelFor i in 0 ..< nBodies.int:
            captures: {bodies, dt, acc_ptr}
            bodies[i].vel += 0.5 * dt * (bodies[i].memory.vector + acc_ptr[i])
            bodies[i].memory.vector = acc_ptr[i]