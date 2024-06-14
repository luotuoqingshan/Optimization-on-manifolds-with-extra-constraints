# for testing
graphs = ["G$i" for i=1:9] 

function RALM_test(graph, seed)::String
    return "ulimit -d $((16 * 1024 * 1024)); "*
        "matlab -singleCompThread -batch"*
        " \"Boss_1_BC(graph='G1', seed=0, solver='RALM');"*
        "Boss_1_BC(graph='$graph', seed=$seed, solver='RALM');\""
end


function Q_LQH_test(graph, seed)::String 
    return "ulimit -d $((16 * 1024 * 1024)); "*
        "matlab -singleCompThread -batch"*
        " \"Boss_1_BC(graph='G1', seed=0, solver='Q_LQH');"*
        "Boss_1_BC(graph='$graph', seed=$seed, solver='Q_LQH');\""
end


function Q_LSE_test(graph, seed)::String 
    return "ulimit -d $((16 * 1024 * 1024)); "*
        "matlab -singleCompThread -batch"*
        " \"Boss_1_BC(graph='G1', seed=0, solver='Q_LSE');"*
        "Boss_1_BC(graph='$graph', seed=$seed, solver='Q_LSE');\""
end

seed = 0
open("test_minimumbisection.txt", "w") do io
    for graph in graphs
        println(io, RALM_test(graph, seed)) 
    end
    for graph in graphs
        println(io, Q_LQH_test(graph, seed)) 
    end
    for graph in graphs
        println(io, Q_LSE_test(graph, seed)) 
    end
end
