using NetworkDynamics

using NetworkDynamics: SpinTask

function busysleep(s)
    starttime = time()
    a = Ref(rand())
    while time() < starttime + s
        a[] = sin(a[])*cos(a[])+rand()
    end
    s
end

t1 = SpinTask("sleep 1") do io
    busysleep(5)
    println(io, "slept 1 second")
    return "from sleep1"
end
t2 = SpinTask("sleep 2") do io
    busysleep(3)
    println(io, "slept 2 seconds")
    return "from sleep2"
end
t3 = SpinTask("sleep 3") do io
    busysleep(3.5)
    # println(io, "slept 3 second")
    return "from sleep3"
end
t4 = SpinTask("sleep 4") do io
    busysleep(4)
    println(io, "slept 4 seconds")
    return "from sleep4"
end
t5 = SpinTask("error task") do io
    busysleep(2.5)
    println(io, "this will error")
    error("intentional error")
end

tasks = [t1, t2, t3, t5, t4]
tasks_noerr = [t1, t2, t3, t4]

NetworkDynamics.run_plain(tasks)
NetworkDynamics.run_plain(tasks_noerr)

NetworkDynamics.run_fancy(tasks)

NetworkDynamics.run_fancy(tasks_noerr)

tasks = [SpinTask((io)->busysleep(rand(1:0.1:2)), "Task $i") for i in 1:30]
NetworkDynamics.run_fancy(tasks)
