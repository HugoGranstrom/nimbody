import tables, math

type
    Vector3* = object
        x*, y*, z*: float64
    Body* = object
        pos*, vel*: Vector3
        mass*: float64
        memory*: tuple[scalar: float64, vector: Vector3]
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