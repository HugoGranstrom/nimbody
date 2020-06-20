import math, tables, random, std/monotimes, times

type
    Vector3* = object
        x*, y*, z*: float64
    Body* = object
        pos*, vel*: Vector3
        mass*: float64
        memory*: tuple[scalar: float64, vector: Vector3]
        #memory*: Table[string, float64]
        #vectorMemory*: Table[string, Vector3]
    OcTree* = ref object
        children*: seq[OcTree]
        bodies*: seq[ptr Body]
        side*: float64
        corner*: Vector3 # lower, left, out corner. if box of side s is centered on (0,0,0) this corner is at (-s/2, -s/2, -s/2)
        centerOfMass*: Vector3
        mass*: float64
        isLeaf*: bool
    RawOcTree* = object
        children*: seq[RawOcTree]
        bodies*: seq[ptr Body]
        side*: float64
        corner*: Vector3 # lower, left, out corner. if box of side s is centered on (0,0,0) this corner is at (-s/2, -s/2, -s/2)
        centerOfMass*: Vector3
        mass*: float64
        isLeaf*: bool
    PersistantMemory*[T] = ref object
        memory*: Table[string, T]
    CalcAccProc*[T, U, V] = proc(timestep: int64, bodies: ptr UncheckedArray[Body], nBodies: int64, floatMem: PersistantMemory[float64], customMem1: PersistantMemory[T], customMem2: PersistantMemory[U], customMem3: PersistantMemory[V]): seq[Vector3]
    IntegratorProc*[T, U, V] = proc(timestep: int64, bodies: ptr UncheckedArray[Body], nBodies: int64, calcAcc: CalcAccProc[T, U, V], dt: float64, floatMem: PersistantMemory[float64], customMem1: PersistantMemory[T], customMem2: PersistantMemory[U], customMem3: PersistantMemory[V])


proc `+`*(v1, v2: Vector3): Vector3 {.inline.} =
    result.x = v1.x + v2.x
    result.y = v1.y + v2.y
    result.z = v1.z + v2.z

proc `+=`*(v1: var Vector3, v2: Vector3) {.inline.} =
    v1.x += v2.x
    v1.y += v2.y
    v1.z += v2.z

proc `-`*(v1, v2: Vector3): Vector3 {.inline.} =
    result.x = v1.x - v2.x
    result.y = v1.y - v2.y
    result.z = v1.z - v2.z

proc `*`*(d: float64, v1: Vector3): Vector3 {.inline.} =
    result.x = d * v1.x
    result.y = d * v1.y
    result.z = d * v1.z

proc `*`*(v1: Vector3, d: float64): Vector3 {.inline.} =
    result.x = d * v1.x
    result.y = d * v1.y
    result.z = d * v1.z

proc `/`*(v1: Vector3, d: float64): Vector3 {.inline.} =
    result.x = v1.x / d
    result.y = v1.y / d
    result.z = v1.z / d

proc mag*(v1: Vector3): float64 {.inline.} =
    result = sqrt(v1.x * v1.x + v1.y * v1.y + v1.z * v1.z)

proc dot*(v1, v2: Vector3): float64 {.inline.} =
    result = v1.x * v2.x + v1.y * v2.y + v1.z * v2.z


proc len*(tree: OcTree): int {.inline.} =
    tree.children.len

proc `$`*(tree: OcTree): string {.inline.} =
    $tree[]

proc isIn*(pos: Vector3, tree: OcTree): bool {.inline.} =
    let s = tree.side
    let corner = tree.corner
    if corner.x < pos.x and pos.x < corner.x + s and corner.y < pos.y and pos.y < corner.y + s and corner.z < pos.z and pos.z < corner.z + s:
        return true
    else:
        return false

proc isIn*(pos: Vector3, corner: Vector3, s: float64): bool {.inline.} =
    if corner.x < pos.x and pos.x < corner.x + s and corner.y < pos.y and pos.y < corner.y + s and corner.z < pos.z and pos.z < corner.z + s:
        return true
    else:
        return false

proc childCorners*(tree: OcTree): array[8, Vector3] {.inline.} =
    let ds = tree.side/2
    result[0] = tree.corner
    result[1] = tree.corner + Vector3(y: ds)
    result[2] = tree.corner + Vector3(x: ds)
    result[3] = tree.corner + Vector3(x: ds, y: ds)
    result[4] = tree.corner + Vector3(z: ds)
    result[5] = tree.corner + Vector3(y: ds, z: ds)
    result[6] = tree.corner + Vector3(x: ds, z: ds)
    result[7] = tree.corner + Vector3(x: ds, y: ds, z: ds)


proc isIn*(pos: Vector3, tree: ptr RawOcTree): bool {.inline.} =
    let s = tree[].side
    let corner = tree[].corner
    if corner.x < pos.x and pos.x < corner.x + s and corner.y < pos.y and pos.y < corner.y + s and corner.z < pos.z and pos.z < corner.z + s:
        return true
    else:
        return false

proc childCorners*(tree: ptr RawOcTree): array[8, Vector3] {.inline.} =
    let ds = tree[].side/2
    result[0] = tree[].corner
    result[1] = tree[].corner + Vector3(y: ds)
    result[2] = tree[].corner + Vector3(x: ds)
    result[3] = tree[].corner + Vector3(x: ds, y: ds)
    result[4] = tree[].corner + Vector3(z: ds)
    result[5] = tree[].corner + Vector3(y: ds, z: ds)
    result[6] = tree[].corner + Vector3(x: ds, z: ds)
    result[7] = tree[].corner + Vector3(x: ds, y: ds, z: ds)

proc `[]=`*[T](mem: PersistantMemory[T], key: string, val: T) {.inline.} =
    mem.memory[key] = val

proc `[]`*[T](mem: PersistantMemory[T], key: string): T {.inline.} =
    mem.memory[key]

proc `$`*[T](mem: PersistantMemory[T]): string =
    $mem.memory
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

proc calcGravityAcc*[T, U, V](timestep: int64, bodies: ptr UncheckedArray[Body], nBodies: int64, floatMem: PersistantMemory[float64], customMem1: PersistantMemory[T], customMem2: PersistantMemory[U], customMem3: PersistantMemory[V]): seq[Vector3] =
    const G = 1e-3
    result = newSeq[Vector3](nBodies)
    for i in 0 ..< nBodies:
        let body1 = bodies[i]
        var acc = Vector3()
        for j in 0 ..< nBodies:
            if i == j:
                continue
            let body2 = bodies[j]
            let r_vec = body2.pos - body1.pos
            let r = r_vec.mag
            let factor = G * body2.mass / (r * r * r)
            acc += factor * r_vec
        result[i] = acc

proc recursiveBarnesHut(tree: OcTree) {.inline.} =
    let nBodies = tree.bodies.len
    if nBodies == 1:
        tree.isLeaf = true
        tree.centerOfMass = tree.bodies[0][].pos
        tree.mass = tree.bodies[0][].mass
        return
    var bodies_array: array[8, seq[ptr Body]]
    let corners = childCorners(tree)
    var mass: float64 = 0.0
    var com = Vector3()
    for i in 0 ..< nBodies:
        for j in 0 ..< corners.len:
            if tree.bodies[i][].pos.isIn(corners[j], tree.side/2):
                bodies_array[j].add(tree.bodies[i])
                break
        mass += tree.bodies[i][].mass
        com += tree.bodies[i][].mass * tree.bodies[i][].pos
    com = com / mass
    tree.mass = mass
    tree.centerOfMass = com
    for i in 0 ..< corners.len:
        if bodies_array[i].len != 0:
            let childTree = OcTree(side: tree.side/2, corner: corners[i], bodies: bodies_array[i])
            tree.children.add(childTree)
            recursiveBarnesHut(childTree)

proc calcAccTreeBody*(tree: Octree, body: Body, G: float64): Vector3 {.inline.} =
    let r_vec = tree.centerOfMass - body.pos
    let r = mag(r_vec)
    let factor = G * tree.mass / (r * r * r)
    result = factor * r_vec

proc calcAccTreeBody*(tree: Octree, body: Body, G: float64, r_vec: Vector3, r: float64): Vector3 {.inline.} =
    let factor = G * tree.mass / (r * r * r)
    result = factor * r_vec


proc barnesHutAccOnBody*(tree: Octree, body: Body, G, theta: float64): Vector3 =
    if tree.isLeaf and isIn(body.pos, tree):
        return Vector3()
    elif tree.isLeaf:
        return calcAccTreeBody(tree, body, G)
    elif isIn(body.pos, tree):
        for i in 0 ..< tree.children.len:
            result += barnesHutAccOnBody(tree.children[i], body, G, theta)
        return
    let r_vec = tree.centerOfMass - body.pos
    let r = mag(r_vec)
    if tree.side / r < theta:
        return calcAccTreeBody(tree, body, G, r_vec, r)
    for i in 0 ..< tree.children.len:
        result += barnesHutAccOnBody(tree.children[i], body, G, theta)
    


proc barnesHutAcc*[T, U, V](timestep: int64, bodies: ptr UncheckedArray[Body], nBodies: int64, floatMem: PersistantMemory[float64], customMem1: PersistantMemory[T], customMem2: PersistantMemory[U], customMem3: PersistantMemory[V]): seq[Vector3] =
    const max_side = 1000
    const G = 1e-3
    const theta = 0.8
    result = newSeq[Vector3](nBodies)
    let root = OcTree(side: max_side, isLeaf: false, corner: Vector3(x: -max_side/2, y: -max_side/2, z: -max_side/2), children: newSeqOfCap[OcTree](8), bodies: newSeqOfCap[ptr Body](nBodies))
    #let corners = childCorners(root)
    #for i in 0 ..< corners.len:
    #    root.children[i] = OcTree(side: root.side/2, corner: corners[i], children: newSeqOfCap[OcTree](8))
    for i in 0 ..< nBodies:
        if bodies[i].pos.isIn(root):
            root.bodies.add(bodies[i].unsafeAddr)
        else:
            when false:
                echo bodies[i].pos, " is out of bounds!"
            else:
                discard
    recursiveBarnesHut(root) # construct tree
    #echo root
    # calculate forces on all bodies:
    for i in 0 ..< nBodies:
        result[i] = barnesHutAccOnBody(root, bodies[i], G, theta)


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

when compileOption("threads"):
    import weave
    # define all built-in integrators but parallel variants in a when statement


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


    proc calcGravityAccParallel*[T, U, V](timestep: int64, bodies: ptr UncheckedArray[Body], nBodies: int64, floatMem: PersistantMemory[float64], customMem1: PersistantMemory[T], customMem2: PersistantMemory[U], customMem3: PersistantMemory[V]): seq[Vector3] =
        result = newSeq[Vector3](nBodies)
        let result_ptr = cast[ptr UncheckedArray[Vector3]](result[0].unsafeAddr)
        const G = 1e-3
        parallelFor i in 0 ..< nBodies.int:
            captures: {G, result_ptr, bodies, nBodies}
            let body1 = bodies[i]
            var acc = Vector3()
            for j in 0 ..< nBodies:
                if i == j:
                    continue
                let body2 = bodies[j]
                let r_vec = body2.pos - body1.pos
                let r = r_vec.mag
                let factor = G * body2.mass / (r * r * r)
                acc += factor * r_vec
            result_ptr[i] = acc
        syncRoot(Weave)
    
    proc velocityVerletStepParallel*[T, U, V](timestep: int64, bodies: ptr UncheckedArray[Body], nBodies: int64, calcAcc: CalcAccProc[T, U, V], dt: float64, floatMem: PersistantMemory[float64], customMem1: PersistantMemory[T], customMem2: PersistantMemory[U], customMem3: PersistantMemory[V]) =
        let dt = dt
        if timestep == 0:
            let accFirst = calcAcc(timestep, bodies, nBodies, floatMem, customMem1, customMem2, customMem3)
            let acc_ptr = cast[ptr UncheckedArray[Vector3]](accFirst[0].unsafeAddr)
            syncScope():
                parallelFor i in 0 ..< nBodies.int:
                    captures: {bodies, acc_ptr}
                    bodies[i].memory.vector = acc_ptr[i]
            #syncRoot(Weave)
        #syncRoot(Weave)
        syncScope():
            parallelFor i in 0 ..< nBodies.int:
                captures: {bodies, dt}
                bodies[i].pos += bodies[i].vel * dt + 0.5 * bodies[i].memory.vector * dt * dt
        #syncRoot(Weave)
        let acc = calcAcc(timestep, bodies, nBodies, floatMem, customMem1, customMem2, customMem3)
        let acc_ptr = cast[ptr UncheckedArray[Vector3]](acc[0].unsafeAddr)
        syncScope():
            parallelFor i in 0 ..< nBodies.int:
                captures: {bodies, dt, acc_ptr}
                bodies[i].vel += 0.5 * dt * (bodies[i].memory.vector + acc_ptr[i])
                bodies[i].memory.vector = acc_ptr[i]
        #syncRoot(Weave)
    

    proc recursiveBarnesHutParallel(tree: ptr RawOcTree) {.inline.} =
        let nBodies = tree[].bodies.len
        if nBodies == 1:
            tree[].isLeaf = true
            tree[].centerOfMass = tree[].bodies[0][].pos
            tree[].mass = tree[].bodies[0][].mass
            return
        var bodies_array: array[8, seq[ptr Body]]
        for i in 0 ..< 8:
            bodies_array[i] = newSeqOfCap[ptr Body](nBodies)
        let corners = childCorners(tree)
        var mass: float64 = 0.0
        var com = Vector3()
        for i in 0 ..< nBodies:
            for j in 0 ..< corners.len:
                if tree[].bodies[i][].pos.isIn(corners[j], tree[].side/2):
                    bodies_array[j].add(tree[].bodies[i])
                    break
            mass += tree[].bodies[i][].mass
            com += tree[].bodies[i][].mass * tree[].bodies[i][].pos
        com = com / mass
        tree[].mass = mass
        tree[].centerOfMass = com
        for i in 0 ..< corners.len:
            if bodies_array[i].len != 0:
                let childTree = RawOcTree(side: tree[].side/2, corner: corners[i], bodies: bodies_array[i])
                tree[].children.add(childTree)
        for i in 0 ..< tree[].children.len:
            #captures: {tree}
            recursiveBarnesHutParallel(tree[].children[i].unsafeAddr)
        

    proc calcAccTreeBodyParallel*(tree: ptr RawOcTree, body: Body, G: float64): Vector3 {.inline.} =
        let r_vec = tree[].centerOfMass - body.pos
        let r = mag(r_vec)
        let factor = G * tree[].mass / (r * r * r)
        result = factor * r_vec

    proc calcAccTreeBodyParallel*(tree: ptr RawOcTree, body: Body, G: float64, r_vec: Vector3, r: float64): Vector3 {.inline.} =
        let factor = G * tree[].mass / (r * r * r)
        result = factor * r_vec


    proc barnesHutAccOnBodyParallel*(tree: ptr RawOcTree, body: Body, G, theta: float64): Vector3 =
        if tree[].isLeaf and isIn(body.pos, tree):
            return Vector3()
        elif tree[].isLeaf:
            return calcAccTreeBodyParallel(tree, body, G)
        elif isIn(body.pos, tree):
            for i in 0 ..< tree[].children.len:
                result += barnesHutAccOnBodyParallel(tree[].children[i].unsafeAddr, body, G, theta)
            return
        let r_vec = tree[].centerOfMass - body.pos
        let r = mag(r_vec)
        if tree[].side / r < theta:
            return calcAccTreeBodyParallel(tree, body, G, r_vec, r)
        for i in 0 ..< tree.children.len:
            result += barnesHutAccOnBodyParallel(tree[].children[i].unsafeAddr, body, G, theta)
        
    proc barnesHutAccParallel*[T, U, V](timestep: int64, bodies: ptr UncheckedArray[Body], nBodies: int64, floatMem: PersistantMemory[float64], customMem1: PersistantMemory[T], customMem2: PersistantMemory[U], customMem3: PersistantMemory[V]): seq[Vector3] =
        const max_side = 1000
        const G = 1e-3
        const theta = 0.8
        result = newSeq[Vector3](nBodies)
        var root = RawOcTree(side: max_side, isLeaf: false, corner: Vector3(x: -max_side/2, y: -max_side/2, z: -max_side/2), children: newSeqOfCap[RawOcTree](8), bodies: newSeqOfCap[ptr Body](nBodies))
        let root_ptr = root.unsafeAddr
        #let corners = childCorners(root)
        #for i in 0 ..< corners.len:
        #    root.children[i] = OcTree(side: root.side/2, corner: corners[i], children: newSeqOfCap[OcTree](8))
        for i in 0 ..< nBodies:
            if bodies[i].pos.isIn(root.unsafeAddr):
                root.bodies.add(bodies[i].unsafeAddr)
            else:
                when false:
                    echo bodies[i].pos, " is out of bounds!"
                else:
                    discard
        syncScope():
            recursiveBarnesHutParallel(root_ptr) # construct tree
        #echo root
        # calculate forces on all bodies:
        let result_ptr = cast[ptr UncheckedArray[Vector3]](result[0].unsafeaddr)
        syncScope():
            parallelFor i in 0 ..< nBodies.int:
                captures: {result_ptr, root_ptr, bodies, G, theta}
                result_ptr[i] = barnesHutAccOnBodyParallel(root_ptr, bodies[i], G, theta)
else:
    echo "--threads option isn't set so no parallel procs are available"


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
    customMem2Internal["charge-val"] = 2.34
    var bodies = @initBodies # make a copy and convert to seq
    let nBodies = bodies.len
    let bodies_ptr = cast[ptr UncheckedArray[Body]](bodies[0].unsafeaddr) # make an uncheckedarray out of it
    let n_whole_steps = abs(floor(simul_time / dt)).int64 # do n full steps and then a correction-step
    when compileOption("threads"):
        init(Weave)
    let startTime = getMonoTime()
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
    #echo last_dt
    integrator(-1, bodies_ptr, bodies.len, calcAcc, last_dt, floatMem, customMem1Internal, customMem2Internal, customMem3Internal)
    # timestep = -1 on the last iteration
    #kin = kineticEnergyParallel(bodies_ptr, nBodies)
    #pot = gravityPotentialEnergyParallel(bodies_ptr, nBodies, 1e-3)
    #echo "After: Kinetic = ", kin, " Potential = ", pot, " Total = ", kin+pot
    #echo "diff: ", tot - (kin + pot)
    let endTime = getMonoTime()
    let elapsed = inMilliSeconds(endTime - startTime)
    echo "Elapsed time: ", elapsed, " ms"
    when compileOption("threads"):
        exit(Weave)
        
    #echo "Energy diff: ", tot - (kin+pot)
    return (bodies, floatMem, customMem1Internal, customMem2Internal, customMem3Internal)
    
when isMainModule:
    var a = Body(mass: 1.0)
    let b = PersistantMemory[float64]()
    let N = 1_000_000
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
    let (bodiesFinal1, floatM1, M11, M21, M31) = nbodySimulation[float64, float64, int64](initBodies = bodies, calcAcc = barnesHutAccParallel, integrator = velocityVerletStep, simul_time = 1e0, dt = 1, customMem1 = b)
    let (bodiesFinal2, floatM2, M12, M22, M32) = nbodySimulation[float64, float64, int64](initBodies = bodies, calcAcc = barnesHutAcc, integrator = velocityVerletStep, simul_time = 1e0, dt = 1, customMem1 = b)
    #let (bodiesFinal3, floatM3, M13, M23, M33) = nbodySimulation[float64, float64, int64](initBodies = bodies, calcAcc = barnesHutAcc, integrator = velocityVerletStepParallel, simul_time = 1e0, dt = 1e0, customMem1 = b)
    let b1 = bodiesFinal1[0]
    let b2 = bodiesFinal1[0]
    #echo bodies[0]
    #echo b1
    #echo b2
    echo "b1 == b2: ", b1 == b2
    #echo "b1 == b0: ", b1 == bodies[0]
    
    #echo M2[]
    #echo M2
    #echo bodies
    #echo sizeof(c[0])