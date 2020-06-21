import weave
import types
# Naive Methods

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




# Serial Barnes Hut

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
    const theta = 0.5
    result = newSeq[Vector3](nBodies)
    let root = OcTree(side: max_side, isLeaf: false, corner: Vector3(x: -max_side/2, y: -max_side/2, z: -max_side/2), children: newSeqOfCap[OcTree](8), bodies: newSeqOfCap[ptr Body](nBodies))
    for i in 0 ..< nBodies:
        if bodies[i].pos.isIn(root):
            root.bodies.add(bodies[i].unsafeAddr)
        else:
            when false:
                echo bodies[i].pos, " is out of bounds!"
            else:
                discard
    recursiveBarnesHut(root) # construct tree
    # calculate forces on all bodies:
    for i in 0 ..< nBodies:
        result[i] = barnesHutAccOnBody(root, bodies[i], G, theta)




# Parallel Barnes Hut

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
    const theta = 0.5
    result = newSeq[Vector3](nBodies)
    var root = RawOcTree(side: max_side, isLeaf: false, corner: Vector3(x: -max_side/2, y: -max_side/2, z: -max_side/2), children: newSeqOfCap[RawOcTree](8), bodies: newSeqOfCap[ptr Body](nBodies))
    let root_ptr = root.unsafeAddr
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
    # calculate forces on all bodies:
    let result_ptr = cast[ptr UncheckedArray[Vector3]](result[0].unsafeaddr)
    syncScope():
        parallelFor i in 0 ..< nBodies.int:
            captures: {result_ptr, root_ptr, bodies, G, theta}
            result_ptr[i] = barnesHutAccOnBodyParallel(root_ptr, bodies[i], G, theta)