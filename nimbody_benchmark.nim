import std/[monotimes, times, strformat, random]
import weave
import nimbody
when compileOption("threads"):
    proc benchmark() =
        var N = 1000
        var bodies1000 = newSeq[Body](N)
        for i in 0 ..< N:
            let x = (rand(1.0) - 0.5) * 20
            let y = (rand(1.0) - 0.5) * 20
            let z = (rand(1.0) - 0.5) * 20
            let pos = Vector3(x: x, y: y, z: z)
            let vel = Vector3(x: -y, y: x, z: z) / 1000
            let mass = rand(1e-2)
            bodies1000[i] = Body(pos: pos, vel: vel, mass: mass)
        N = 10_000
        var bodies10_000 = newSeq[Body](N)
        for i in 0 ..< N:
            let x = (rand(1.0) - 0.5) * 20
            let y = (rand(1.0) - 0.5) * 20
            let z = (rand(1.0) - 0.5) * 20
            let pos = Vector3(x: x, y: y, z: z)
            let vel = Vector3(x: -y, y: x, z: z) / 1000
            let mass = rand(1e-2)
            bodies10_000[i] = Body(pos: pos, vel: vel, mass: mass)
        N = 100_000
        var bodies100_000 = newSeq[Body](N)
        for i in 0 ..< N:
            let x = (rand(1.0) - 0.5) * 20
            let y = (rand(1.0) - 0.5) * 20
            let z = (rand(1.0) - 0.5) * 20
            let pos = Vector3(x: x, y: y, z: z)
            let vel = Vector3(x: -y, y: x, z: z) / 1000
            let mass = rand(1e-2)
            bodies100_000[i] = Body(pos: pos, vel: vel, mass: mass)
        let simul_time = 1.0
        let dt = 0.3
        let n_iters = 5
        init(Weave)
        echo "Benchmarks for 1000 bodies:"
        var startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies1000, calcAcc = barnesHutAccParallel, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        var endTime = getMonoTime()
        var elipsedBHParallel = inMilliseconds(endTime - startTime)
        echo &"BarnesParallel: {elipsedBHParallel} ms ({elipsedBHParallel.float / n_iters.float} ms/iter) ({elipsedBHParallel.float / elipsedBHParallel.float}x BarnesParallel)"

        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies1000, calcAcc = barnesHutAcc, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedBHSerial= inMilliseconds(endTime - startTime)
        echo &"BarnesSerial: {elipsedBHSerial} ms ({elipsedBHSerial.float / n_iters.float} ms/iter) ({elipsedBHSerial.float / elipsedBHParallel.float}x BarnesParallel)"

        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies1000, calcAcc = calcGravityAccParallel, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedNaiveParallel = inMilliseconds(endTime - startTime)
        echo &"NaiveParallel: {elipsedNaiveParallel} ms ({elipsedNaiveParallel.float / n_iters.float} ms/iter) ({elipsedNaiveParallel.float / elipsedBHParallel.float}x BarnesParallel)"
        
        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies1000, calcAcc = calcGravityAcc, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedNaiveSerial = inMilliseconds(endTime - startTime)
        echo &"NaiveSerial: {elipsedNaiveSerial} ms ({elipsedNaiveSerial.float / n_iters.float} ms/iter) ({elipsedNaiveSerial.float / elipsedBHParallel.float}x BarnesParallel)"
        

        echo "\nBenchmarks for 10_000 bodies:"
        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies10000, calcAcc = barnesHutAccParallel, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedBHParallel2 = inMilliseconds(endTime - startTime)
        echo &"BarnesParallel: {elipsedBHParallel2} ms ({elipsedBHParallel2.float / n_iters.float} ms/iter) ({elipsedBHParallel2.float / elipsedBHParallel2.float}x BarnesParallel)"

        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies10000, calcAcc = barnesHutAcc, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedBHSerial2 = inMilliseconds(endTime - startTime)
        echo &"BarnesSerial: {elipsedBHSerial2} ms ({elipsedBHSerial2.float / n_iters.float} ms/iter) ({elipsedBHSerial2.float / elipsedBHParallel2.float}x BarnesParallel)"

        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies10000, calcAcc = calcGravityAccParallel, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedNaiveParallel2 = inMilliseconds(endTime - startTime)
        echo &"NaiveParallel: {elipsedNaiveParallel2} ms ({elipsedNaiveParallel2.float / n_iters.float} ms/iter) ({elipsedNaiveParallel2.float / elipsedBHParallel2.float}x BarnesParallel)"
        
        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies10000, calcAcc = calcGravityAcc, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedNaiveSerial2 = inMilliseconds(endTime - startTime)
        echo &"NaiveSerial: {elipsedNaiveSerial2} ms ({elipsedNaiveSerial2.float / n_iters.float} ms/iter) ({elipsedNaiveSerial2.float / elipsedBHParallel2.float}x BarnesParallel)"
        

        echo "\nBenchmarks for 100_000 bodies:"
        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies100_000, calcAcc = barnesHutAccParallel, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedBHParallel3 = inMilliseconds(endTime - startTime)
        echo &"BarnesParallel: {elipsedBHParallel3} ms ({elipsedBHParallel3.float / n_iters.float} ms/iter) ({elipsedBHParallel3.float / elipsedBHParallel3.float}x BarnesParallel)"

        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies100_000, calcAcc = barnesHutAcc, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedBHSerial3 = inMilliseconds(endTime - startTime)
        echo &"BarnesSerial: {elipsedBHSerial3} ms ({elipsedBHSerial3.float / n_iters.float} ms/iter) ({elipsedBHSerial3.float / elipsedBHParallel3.float}x BarnesParallel)"

        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies100_000, calcAcc = calcGravityAccParallel, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedNaiveParallel3 = inMilliseconds(endTime - startTime)
        echo &"NaiveParallel: {elipsedNaiveParallel3} ms ({elipsedNaiveParallel3.float / n_iters.float} ms/iter) ({elipsedNaiveParallel3.float / elipsedBHParallel3.float}x BarnesParallel)"
        
        startTime = getMonoTime()
        discard nbodySimulation[float64, float64, float64](initBodies = bodies100_000, calcAcc = calcGravityAcc, integrator = velocityVerletStep, simul_time = simul_time, dt = dt)
        endTime = getMonoTime()
        let elipsedNaiveSerial3 = inMilliseconds(endTime - startTime)
        echo &"NaiveSerial: {elipsedNaiveSerial3} ms ({elipsedNaiveSerial3.float / n_iters.float} ms/iter) ({elipsedNaiveSerial3.float / elipsedBHParallel3.float}x BarnesParallel)"
        
        
        exit(Weave)

    benchmark()