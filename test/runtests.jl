using Test
using CoNCMOR

module mpart3
using StaticArrays
using SparseArrays
using Test
# using PlotlyJS
using FinEtools
using CoNCMOR: CoNCData, nclusters, nfuncspercluster
using CoNCMOR: nbasisfunctions, transfmatrix, LegendreBasis
colors = ["rgb(164, 194, 244)", "rgb(194, 194, 144)", "rgb(194, 144, 244)", "rgb(164, 244, 144)", "rgb(164, 194, 244)", "rgb(255, 217, 102)", "rgb(234, 153, 153)", "rgb(142, 124, 195)"]
function test()
    points = StaticArrays.SArray{Tuple{2},Float64,1,2}[[0.6501312860170676, 0.889978628556229], [0.19162369745368268, 0.5000498373796463], [0.9513777178320113, 0.9313997945413548], [0.6715708752191747, 0.0031159439825900748], [0.908983887358429, 0.33697103964932174], [0.8630730459099334, 0.8172618688642774], [0.21191150611345178, 0.6279042354095665], [0.0017150255179243512, 0.7401108604352873], [0.045937411259845184, 0.3867144891772889], [0.07794203801162314, 0.45466525770422384], [0.11103905573968631, 0.7213589211503646], [0.522501483190954, 0.9505399568028736], [0.14011031960881337, 0.5976293299997106], [0.15824295891177775, 0.4538838028009915], [0.22532364404660643, 0.7490903677067646], [0.16495650937348727, 0.4095197672868107], [0.2929326909660668, 0.2627344309574493], [0.7749630590488474, 0.010281406373283675], [0.23137166232629092, 0.1872554576301051], [0.48014323786400226, 0.8846143614245845], [0.26528544579757507, 
    0.6924068278632542], [0.2839884809072999, 0.09439141875666568], [0.942287809602099, 0.4005321413670275], [0.026409937646382442, 0.6845508696068237], [0.7752785557256849, 0.3285662200502397], [0.09380966727246509, 0.15933031973012612], [0.2813301059379758, 0.15216808040457552], [0.4780694864591679, 0.6063568972874278], [0.6982811071004462, 0.5308171330635794], [0.9609668367122937, 0.09387988328210572]]  
    X = fill(0.0, length(points), 2)
    for i in 1:length(points)
        X[i, :] .= points[i]
    end        
    ppartitioning = pointpartitioning(X, 2)
    partitionnumbers = unique(ppartitioning)

    # data = PlotlyBase.AbstractTrace[]
    # for gp in partitionnumbers
    #     trace1 = scatter(; 
    #         x=[points[i][1] for i in 1:length(points) if ppartitioning[i] == gp], 
    #         y=[points[i][2] for i in 1:length(points) if ppartitioning[i] == gp] ,
    #       mode="markers",
    #       marker=attr(color=colors[gp], size=12,
    #           line=attr(color="white", width=0.5))
    #       )
    #     push!(data, trace1)
    # end
    # layout = Layout(;title="Point partitioning",
    #     xaxis=attr(title="x", zeroline=false),
    #     yaxis=attr(title="y", zeroline=false))

    # pl = plot(data, layout)
    # display(pl)

    mor = CoNCData(X, ppartitioning)
    @test nclusters(mor) == length(partitionnumbers)
    geom = NodalField(X) 
    u = NodalField(X) 
    numberdofs!(u)
    Phi  = transfmatrix(mor, LegendreBasis, 1, u)
    @test nfuncspercluster(mor) == [1, 1]
    @test nbasisfunctions(mor) == 2
    
    I, J, V = findnz(Phi)
    @test (I, J, V) == ([1, 5, 7, 9, 11, 23, 35, 45, 49, 55, 57, 59, 2, 6, 8, 10, 12, 24, 36, 46, 50, 56, 58, 60, 3, 13, 15, 17, 19, 21, 25, 27, 29, 31, 33, 37, 39, 41, 43, 47, 51, 53, 4, 14, 16, 18, 20, 22, 26, 28, 30, 32, 34, 38, 40, 42, 44, 48, 52, 54], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
    4, 4, 4, 4, 4, 4, 4], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    nothing
end
test()
end



module mpart4
using StaticArrays
using SparseArrays
using Test
# using PlotlyJS
using FinEtools
using CoNCMOR: CoNCData, nclusters, nfuncspercluster
using CoNCMOR: nbasisfunctions, transfmatrix, LegendreBasis
colors = ["rgb(164, 194, 244)", "rgb(194, 194, 144)", "rgb(194, 144, 244)", "rgb(164, 244, 144)", "rgb(164, 194, 244)", "rgb(255, 217, 102)", "rgb(234, 153, 153)", "rgb(142, 124, 195)"]
function test()
    points = StaticArrays.SArray{Tuple{2},Float64,1,2}[[0.6501312860170676, 0.889978628556229], [0.19162369745368268, 0.5000498373796463], [0.9513777178320113, 0.9313997945413548], [0.6715708752191747, 0.0031159439825900748], [0.908983887358429, 0.33697103964932174], [0.8630730459099334, 0.8172618688642774], [0.21191150611345178, 0.6279042354095665], [0.0017150255179243512, 0.7401108604352873], [0.045937411259845184, 0.3867144891772889], [0.07794203801162314, 0.45466525770422384], [0.11103905573968631, 0.7213589211503646], [0.522501483190954, 0.9505399568028736], [0.14011031960881337, 0.5976293299997106], [0.15824295891177775, 0.4538838028009915], [0.22532364404660643, 0.7490903677067646], [0.16495650937348727, 0.4095197672868107], [0.2929326909660668, 0.2627344309574493], [0.7749630590488474, 0.010281406373283675], [0.23137166232629092, 0.1872554576301051], [0.48014323786400226, 0.8846143614245845], [0.26528544579757507, 
    0.6924068278632542], [0.2839884809072999, 0.09439141875666568], [0.942287809602099, 0.4005321413670275], [0.026409937646382442, 0.6845508696068237], [0.7752785557256849, 0.3285662200502397], [0.09380966727246509, 0.15933031973012612], [0.2813301059379758, 0.15216808040457552], [0.4780694864591679, 0.6063568972874278], [0.6982811071004462, 0.5308171330635794], [0.9609668367122937, 0.09387988328210572]]  
    X = fill(0.0, length(points), 2)
    for i in 1:length(points)
        X[i, :] .= points[i]
    end        
    ppartitioning = pointpartitioning(X, 2)
    partitionnumbers = unique(ppartitioning)

    # data = PlotlyBase.AbstractTrace[]
    # for gp in partitionnumbers
    #     trace1 = scatter(; 
    #         x=[points[i][1] for i in 1:length(points) if ppartitioning[i] == gp], 
    #         y=[points[i][2] for i in 1:length(points) if ppartitioning[i] == gp] ,
    #       mode="markers",
    #       marker=attr(color=colors[gp], size=12,
    #           line=attr(color="white", width=0.5))
    #       )
    #     push!(data, trace1)
    # end
    # layout = Layout(;title="Point partitioning",
    #     xaxis=attr(title="x", zeroline=false),
    #     yaxis=attr(title="y", zeroline=false))

    # pl = plot(data, layout)
    # display(pl)

    mor = CoNCData(X, ppartitioning)
    @test nclusters(mor) == length(partitionnumbers)
    geom = NodalField(X) 
    u = NodalField(X) 
    numberdofs!(u)
    Phi  = transfmatrix(mor, LegendreBasis, 2, u)
    @test nfuncspercluster(mor) == [3, 3]
    @test nbasisfunctions(mor) == 2*3
    
    # I, J, V = findnz(Phi)
    # @test (I, J, V) == ([1, 5, 7, 9, 11, 23, 35, 45, 49, 55, 57, 59, 2, 6, 8, 10, 12, 24, 36, 46, 50, 56, 58, 60, 3, 13, 15, 17, 19, 21, 25, 27, 29, 31, 33, 37, 39, 41, 43, 47, 51, 53, 4, 14, 16, 18, 20, 22, 26, 28, 30, 32, 34, 38, 40, 42, 44, 48, 52, 54], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 
    # 4, 4, 4, 4, 4, 4, 4], [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    nothing
end
test()
end


module mesh_Q4spheren_1
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData
using Statistics
using Test
function test()
    rex =  2.0; #external radius
    nr = 13;
    npanelgroups = 4
    
    fens, fes = Q4spheren(rex, nr)
    
    vtkexportmesh("sphere.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.VTK.Q4)
    
    femm1  =  FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    partitioning = Metis.partition(g, npanelgroups; alg = :KWAY)
    X = fill(zero(Float64), count(fes), 3)
    bconn = connasarray(fes)
    for i in 1:size(bconn, 1)
        X[i, :] = @views mean(fens.xyz[bconn[i, :], :], dims=1)
    end
    mor = CoNCData(X, partitioning)
    
    @test count(fens) == 169
    @test count(fes) == 147
end
test()
end

#         mor = CoNCData(X, partitioning)
#         P = ElementalField(zeros(size(X,1), 1)) # Pressure field
#         numberdofs!(P) #
#         Phi = transfmatrix(mor, LegendreBasis, nbf1max, P);


module mesh_Q4spheren_2
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData
using Statistics
using Test
function test()
    rex =  2.0; #external radius
    nr = 13;
    npanelgroups = 4
    
    fens, fes = Q4spheren(rex, nr)
    
    vtkexportmesh("sphere.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.VTK.Q4)
    
    femm1  =  FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    partitioning = Metis.partition(g, npanelgroups; alg = :KWAY)
    partitions = unique(partitioning)
    for j in 1:length(partitions)
        sfes = subset(fes, findall(v -> v == partitions[j], partitioning))
        vtkexportmesh("sphere-p$(partitions[j]).vtk", fens, sfes)
    end
    X = fill(zero(Float64), count(fes), 3)
    bconn = connasarray(fes)
    for i in 1:size(bconn, 1)
        X[i, :] = @views mean(fens.xyz[bconn[i, :], :], dims=1)
    end
    function coordinates(list)
        X[list, :]
    end

    mor = CoNCData(coordinates, partitioning)
    
    @test count(fens) == 169
    @test count(fes) == 147
end
test()
end


module mesh_Q4spheren_3
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData
using Statistics
using LinearAlgebra
using Test
function test()
    rex =  2.0; #external radius
    nr = 4;
    npanelgroups = 4
    
    fens, fes = Q4block(2*rex, rex, nr, nr)
    fens.xyz = xyz3(fens)
    
    vtkexportmesh("block.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.VTK.Q4)
    
    femm1  =  FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    partitioning = Metis.partition(g, npanelgroups; alg = :KWAY)
    partitions = unique(partitioning)
    for j in 1:length(partitions)
        sfes = subset(fes, findall(v -> v == partitions[j], partitioning))
        vtkexportmesh("block-p$(partitions[j]).vtk", fens, sfes)
    end
    X = fill(zero(Float64), count(fes), 3)
    bconn = connasarray(fes)
    for i in 1:size(bconn, 1)
        X[i, :] = @views mean(fens.xyz[bconn[i, :], :], dims=1)
    end
    function coordinates(list)
        krondelta(i, k) = i == k ? 1.0 : 0.0
        lX = X[list, :]
        center = mean(lX, dims = 1)
        for j in 1:size(lX, 1)
            lX[j, :] -= center[:]
        end
        It = fill(0.0, 3, 3)
        for j in 1:size(lX, 1)
            r2 = dot(lX[j, :], lX[j, :])
            for i in 1:3
                for k in i:3
                    It[i, k] += krondelta(i, k) * r2 - lX[j, i] * lX[j, k]
                end
            end
        end
        for i in 1:3
            for k in 1:i-1
                It[i, k] = It[k, i]
            end
        end
        @assert It ==  It'
        # @show center
        # @show It
        epsol = eigen(It)
        normal = epsol.vectors[:, 3]
        lX
    end

    mor = CoNCData(coordinates, partitioning)
    
    @test count(fens) == (nr+1)^2
    @test count(fes) == nr^2
end
test()
end


module mesh_Q4spheren_4
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData
using Statistics
using LinearAlgebra
using Test
function test()
    rex =  2.0; #external radius
    nr = 8;
    npanelgroups = 4
    
    fens, fes = Q4spheren(rex, nr)
    
    vtkexportmesh("sphere.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.VTK.Q4)
    
    femm1  =  FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    partitioning = Metis.partition(g, npanelgroups; alg = :KWAY)
    partitions = unique(partitioning)
    @assert length(partitions) == npanelgroups
    for j in 1:length(partitions)
        sfes = subset(fes, findall(v -> v == partitions[j], partitioning))
        vtkexportmesh("sphere-p$(partitions[j]).vtk", fens, sfes)
    end
    X = fill(zero(Float64), count(fes), 3)
    bconn = connasarray(fes)
    for i in 1:size(bconn, 1)
        X[i, :] = @views mean(fens.xyz[bconn[i, :], :], dims=1)
    end
    function coordinates(list)
        krondelta(i, k) = i == k ? 1.0 : 0.0
        lX = X[list, :]
        center = mean(lX, dims = 1)
        for j in 1:size(lX, 1)
            lX[j, :] -= center[:]
        end
        It = fill(0.0, 3, 3)
        for j in 1:size(lX, 1)
            r2 = dot(lX[j, :], lX[j, :])
            for i in 1:3
                for k in i:3
                    It[i, k] += krondelta(i, k) * r2 - lX[j, i] * lX[j, k]
                end
            end
        end
        for i in 1:3
            for k in 1:i-1
                It[i, k] = It[k, i]
            end
        end
        @assert It ==  It'
        epsol = eigen(It)
        normal = epsol.vectors[:, 3]
        e1 = epsol.vectors[:, 1]
        e2 = epsol.vectors[:, 2]
        lxy = fill(0.0, length(list), 2)
        for j in 1:size(lX, 1)
            lxy[j, 1] = dot(lX[j, :], e1)
            lxy[j, 2] = dot(lX[j, :], e2)
        end
        return lxy
    end

    mor = CoNCData(coordinates, partitioning)
    
    # @test count(fens) == (nr+1)^2
    # @test count(fes) == 147
end
test()
end

module mesh_Q4spheren_5
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData, transfmatrix, LegendreBasis
using Statistics
using LinearAlgebra
using Test
function test()
    rex = 2.0 #external radius
    nr = 45
    npanelgroups = 16
    nbf1max = 4
    function pressure(xy)
        sin(4.0*xy[1]) + sin(2.0*xy[3])
    end

    fens, fes = Q4spheren(rex, nr)

    # vtkexportmesh("sphere.vtk", fes.conn, fens.xyz,  FinEtools.MeshExportModule.VTK.Q4)

    femm1 = FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    partitioning = Metis.partition(g, npanelgroups; alg=:KWAY)
    partitions = unique(partitioning)
    @assert length(partitions) == npanelgroups
    for j in 1:length(partitions)
        sfes = subset(fes, findall(v -> v == partitions[j], partitioning))
        vtkexportmesh("sphere-p$(partitions[j]).vtk", fens, sfes)
    end
    X = fill(zero(Float64), count(fes), 3)
    bconn = connasarray(fes)
    for i in 1:size(bconn, 1)
        X[i, :] = @views mean(fens.xyz[bconn[i, :], :], dims=1)
    end
    function coordinates(list)
        krondelta(i, k) = i == k ? 1.0 : 0.0
        lX = X[list, :]
        center = mean(lX, dims=1)
        for j in 1:size(lX, 1)
            lX[j, :] -= center[:]
        end
        It = fill(0.0, 3, 3)
        for j in 1:size(lX, 1)
            r2 = dot(lX[j, :], lX[j, :])
            for i in 1:3
                for k in i:3
                    It[i, k] += krondelta(i, k) * r2 - lX[j, i] * lX[j, k]
                end
            end
        end
        for i in 1:3
            for k in 1:i-1
                It[i, k] = It[k, i]
            end
        end
        @assert It == It'
        epsol = eigen(It)
        normal = epsol.vectors[:, 3]
        e1 = epsol.vectors[:, 1]
        e2 = epsol.vectors[:, 2]
        lxy = fill(0.0, length(list), 2)
        for j in 1:size(lX, 1)
            lxy[j, 1] = dot(lX[j, :], e1)
            lxy[j, 2] = dot(lX[j, :], e2)
        end
        return lxy
    end

    mor = CoNCData(coordinates, partitioning)
    P = ElementalField(zeros(size(X, 1), 1)) # Pressure field
    numberdofs!(P) #
    Phi = transfmatrix(mor, LegendreBasis, nbf1max, P)

    for j in 1:size(X, 1)
        P.values[j] = pressure(X[j, :])
    end 
    vtkexportmesh("sphere-P.vtk", fens, fes; scalars = [("P", deepcopy(P.values))])

    # qrf = qr(Phi)
    # a = qrf.R \ (Matrix(qrf.Q)' * P.values)
    # Pa = Phi * a
    Pa = Phi * ((Phi' * Phi) \ (Phi' * P.values))
    vtkexportmesh("sphere-Pa.vtk", fens, fes; scalars = [("Pa", deepcopy(Pa))])

    vtkexportmesh("sphere-DeltaP.vtk", fens, fes; scalars = [("DeltaP", deepcopy(Pa - P.values))])

    # @test count(fens) == (nr+1)^2
    # @test count(fes) == 147
end
test()
end


module mesh_Q4multi_1
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData, transfmatrix, LegendreBasis
using Statistics
using LinearAlgebra
using Test

function pressure(xy)
    sin(4.0*xy[1]) + sin(2.0*xy[2])
end

function test()
    A, B = 2.0, 3.0
    nA, nB = 45, 55
    nelgroups = 3
    psize = 100
    
    fens, fes = Q4block(A, B, nA, nB)

    femm1 = FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    epartitioning = Metis.partition(g, nelgroups; alg=:KWAY)
    epartitions = unique(epartitioning)

    for j in eachindex(epartitions)
        sfes = subset(fes, findall(v -> v == epartitions[j], epartitioning))
        vtkexportmesh("rblock1-ep$(epartitions[j]).vtk", fens, sfes)
    end

    gpartitioning = fill(0, count(fens))
    allpartitions = Int[]
    cpartoffset = 0
    for j in eachindex(epartitions)
        # @show "Partition $j"
        sfes = subset(fes, findall(v -> v == epartitions[j], epartitioning))
        femm1 = FEMMBase(IntegDomain(sfes, GaussRule(2, 1)))
        C = connectionmatrix(femm1, count(fens))
        g = Metis.graph(C; check_hermitian=true)
        nl = connectednodes(sfes)
        np = Int(round(length(nl) / psize))
        ppartitioning = Metis.partition(g, np; alg=:KWAY)
        ppartitions = unique(ppartitioning[nl])
        nap = length(allpartitions) + length(ppartitions)
        allpartitions = vcat(allpartitions, ppartitions)
        @assert nap == length(allpartitions)
        for k in eachindex(nl)
            gpartitioning[nl[k]] = ppartitioning[nl[k]] + cpartoffset
        end 
        gpartitions = unique(gpartitioning)
        cpartoffset += length(ppartitions)
    end 
    @assert isempty(findall(v -> v == 0, gpartitioning))
    gpartitions = unique(gpartitioning)


    @test length(unique(gpartitioning)) == 12

    X = fens.xyz
    mor = CoNCData(X, gpartitioning)
    P = NodalField(zeros(size(X, 1), 1)) # Pressure field
    numberdofs!(P) #
    Phi = transfmatrix(mor, LegendreBasis, 5, P)

    for j in 1:size(X, 1)
        P.values[j] = pressure(X[j, :])
    end 
    vtkexportmesh("rblock1-P.vtk", fens, fes; scalars = [("P", deepcopy(P.values))])

    Pa = Phi * ((Phi' * Phi) \ (Phi' * P.values))
    vtkexportmesh("rblock1-Pa.vtk", fens, fes; scalars = [("Pa", deepcopy(Pa))])


    for j in 1:length(gpartitions)
        l = findall(v -> v == gpartitions[j], gpartitioning)
        sfes = FESetP1(reshape(l, length(l), 1))
        vtkexportmesh("rblock1-np$(gpartitions[j]).vtk", fens, sfes)
    end


end
test()
end

module mesh_Q4multi_2
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData, transfmatrix, LegendreBasis
using Statistics
using LinearAlgebra
using Test

function pressure(xy)
    sin(4.0*xy[1]) + sin(2.0*xy[2])
end

function test()
    A, B = 2.0, 3.0
    nA, nB = 21, 23
    nelgroups = 3
    psize = 20
    dontpartition = 0
    
    fens, fes = Q4block(A, B, nA, nB)

    # Here we create the domain partitioning based on the finite elements.
    femm1 = FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    dpartitioning = Metis.partition(g, nelgroups; alg=:KWAY)
    dpartitions = unique(dpartitioning)

    for j in eachindex(dpartitions)
        sfes = subset(fes, findall(v -> v == dpartitions[j], dpartitioning))
        vtkexportmesh("rblock2-domain$(dpartitions[j]).vtk", fens, sfes)
    end

    gpartitioning = fill(0, count(fens))
    allpartitions = Int[]
    cpartoffset = 0
    for j in eachindex(dpartitions)
        sfes = subset(fes, findall(v -> v == dpartitions[j], dpartitioning))
        nl = connectednodes(sfes)
        if j == dontpartition
            ppartitioning = fill(0, count(fens))
            ppartitioning[nl] = 1:length(nl)
        else
            femm1 = FEMMBase(IntegDomain(sfes, GaussRule(2, 1)))
            C = connectionmatrix(femm1, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            np = Int(round(length(nl) / psize))
            # This partitions all the nodes, not just connected by the elements
            # above. 
            ppartitioning = Metis.partition(g, np; alg=:KWAY)
            # Therefore, we need to relabel the partition numbers. All the
            # partitions for nodes not connected together by fes will be set to
            # zero.
            epts = ppartitioning[nl]
            ppartitioning = fill(0, count(fens))
            ppartitioning[nl] = epts
            # Now we have the the current partitions. They are not necessarily
            # numbered consecutively. That needs to be fixed now.
            ppartitions = unique(epts)
            newppartitions = fill(0, maximum(ppartitions))
            newppartitions[ppartitions] = 1:length(ppartitions)
            for k in eachindex(ppartitioning)
                if ppartitioning[k] > 0
                    ppartitioning[k] = newppartitions[ppartitioning[k]]
                end 
            end 
        end 
        # These are now the partitions of the nodes in the current domain.
        ppartitions = unique(ppartitioning[nl])
        nap = length(allpartitions) + length(ppartitions)
        allpartitions = vcat(allpartitions, ppartitions)
        @assert nap == length(allpartitions)
        for k in eachindex(nl)
            gpartitioning[nl[k]] = ppartitioning[nl[k]] + cpartoffset
        end 
        cpartoffset += length(ppartitions)
    end 

    # Now we have the global partitioning of the nodes. Check that all nodes are
    # included (i.e. the partition is not 0)
    @assert isempty(findall(v -> v == 0, gpartitioning))
    # This is how many partitions of the nodes we will have
    gpartitions = unique(gpartitioning)

    X = fens.xyz
    mor = CoNCData(X, gpartitioning)

    P = NodalField(zeros(size(X, 1), 1)) # Pressure field
    numberdofs!(P) #
    Phi = transfmatrix(mor, LegendreBasis, 2, P)

    for j in 1:size(X, 1)
        P.values[j] = pressure(X[j, :])
    end 
    vtkexportmesh("rblock2-P.vtk", fens, fes; scalars = [("P", deepcopy(P.values))])

    Pa = Phi * ((Phi' * Phi) \ (Phi' * P.values))
    vtkexportmesh("rblock2-Pa.vtk", fens, fes; scalars = [("Pa", deepcopy(Pa))])


    for j in 1:length(gpartitions)
        l = findall(v -> v == gpartitions[j], gpartitioning)
        sfes = FESetP1(reshape(l, length(l), 1))
        vtkexportmesh("rblock2-np$(gpartitions[j]).vtk", fens, sfes)
    end


end
test()
end

module mesh_Q4multi_3
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData, transfmatrix, LegendreBasis
using Statistics
using LinearAlgebra
using Test

function pressure(xy)
    sin(4.0*xy[1]) + sin(2.0*xy[2])
end

function test()
    A, B = 2.0, 3.0
    nA, nB = 21, 23
    nelgroups = 3
    psize = 20
    dontpartition = 2
    
    fens, fes = Q4block(A, B, nA, nB)

    # Here we create the domain partitioning based on the finite elements.
    femm1 = FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    dpartitioning = Metis.partition(g, nelgroups; alg=:KWAY)
    dpartitions = unique(dpartitioning)

    for j in eachindex(dpartitions)
        sfes = subset(fes, findall(v -> v == dpartitions[j], dpartitioning))
        vtkexportmesh("rblock3-domain$(dpartitions[j]).vtk", fens, sfes)
    end

    gpartitioning = fill(0, count(fens))
    allpartitions = Int[]
    cpartoffset = 0
    for j in eachindex(dpartitions)
        sfes = subset(fes, findall(v -> v == dpartitions[j], dpartitioning))
        nl = connectednodes(sfes)
        if j == dontpartition
            ppartitioning = fill(0, count(fens))
            ppartitions = 1:length(nl)
            ppartitioning[nl] .= ppartitions
        else
            femm1 = FEMMBase(IntegDomain(sfes, GaussRule(2, 1)))
            C = connectionmatrix(femm1, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            np = Int(round(length(nl) / psize))
            # This partitions all the nodes, not just connected by the elements
            # above. 
            ppartitioning = Metis.partition(g, np; alg=:KWAY)
            # Therefore, we need to relabel the partition numbers. All the
            # partitions for nodes not connected together by fes will be set to
            # zero.
            epts = ppartitioning[nl]
            ppartitioning = fill(0, count(fens))
            ppartitioning[nl] = epts
            # Now we have the the current partitions. They are not necessarily
            # numbered consecutively. That needs to be fixed now.
            ppartitions = unique(epts)
            newppartitions = fill(0, maximum(ppartitions))
            newppartitions[ppartitions] = 1:length(ppartitions)
            for k in eachindex(ppartitioning)
                if ppartitioning[k] > 0
                    ppartitioning[k] = newppartitions[ppartitioning[k]]
                end 
            end 
            ppartitions = unique(ppartitioning[nl])
        end 
        # These are now the partitions of the nodes in the current domain.
        nap = length(allpartitions) + length(ppartitions)
        allpartitions = vcat(allpartitions, ppartitions)
        @assert nap == length(allpartitions)
        for k in eachindex(nl)
            gpartitioning[nl[k]] = ppartitioning[nl[k]] + cpartoffset
        end 
        cpartoffset += length(ppartitions)
    end 

    # Now we have the global partitioning of the nodes. Check that all nodes are
    # included (i.e. the partition is not 0)
    @assert isempty(findall(v -> v == 0, gpartitioning))
    # Some partitions may be missing: they are not necessarily in sequential
    # order. We need to remember them.
    gpartitions = unique(gpartitioning)
    newgpartitions = fill(0, maximum(gpartitions))
    newgpartitions[gpartitions] = 1:length(gpartitions)
    for k in eachindex(gpartitioning)
        if gpartitioning[k] > 0
            gpartitioning[k] = newgpartitions[gpartitioning[k]]
        end 
    end 
    gpartitions = unique(gpartitioning)
    # @show gpartitioning

    # that now we will start using the nodal partitioning.
    for k in eachindex(gpartitioning)
        @assert gpartitioning[k] in gpartitions
    end 

    X = fens.xyz
    mor = CoNCData(X, gpartitioning)

    P = NodalField(zeros(size(X, 1), 1)) # Pressure field
    numberdofs!(P) #
    Phi = transfmatrix(mor, LegendreBasis, 1, P)

    for j in 1:size(X, 1)
        P.values[j] = pressure(X[j, :])
    end 
    vtkexportmesh("rblock3-P.vtk", fens, fes; scalars = [("P", deepcopy(P.values))])

    Pa = Phi * ((Phi' * Phi) \ (Phi' * P.values))
    vtkexportmesh("rblock3-Pa.vtk", fens, fes; scalars = [("Pa", deepcopy(Pa))])


    for j in 1:length(gpartitions)
        l = findall(v -> v == gpartitions[j], gpartitioning)
        sfes = FESetP1(reshape(l, length(l), 1))
        vtkexportmesh("rblock3-np$(gpartitions[j]).vtk", fens, sfes)
    end


end
test()
end


module mesh_Q4multi_4
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData, transfmatrix, LegendreBasis
using Statistics
using LinearAlgebra
using Test

function pressure(xy)
    sin(7.5*xy[1]) + sin(5.3*xy[2])
end

function test()
    A, B = 2.0, 3.0
    nA, nB = 51, 73
    nelgroups = 5
    psize = 20
    dontpartition = 2
    
    fens, fes = Q4block(A, B, nA, nB)
    
    # Here we create the domain partitioning based on the finite elements.
    femm1 = FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    dpartitioning = Metis.partition(g, nelgroups; alg=:KWAY)
    dpartitions = unique(dpartitioning)

    for j in eachindex(dpartitions)
        sfes = subset(fes, findall(v -> v == dpartitions[j], dpartitioning))
        vtkexportmesh("rblock-domain$(dpartitions[j]).vtk", fens, sfes)
    end

    gpartitioning = fill(0, count(fens))
    allpartitions = Int[]
    cpartoffset = 0
    for j in eachindex(dpartitions)
        sfes = subset(fes, findall(v -> v == dpartitions[j], dpartitioning))
        nl = connectednodes(sfes)
        if dpartitions[j] == dontpartition
            ppartitioning = fill(0, count(fens))
            ppartitions = 1:length(nl)
            ppartitioning[nl] .= ppartitions
        else
            femm1 = FEMMBase(IntegDomain(sfes, GaussRule(2, 1)))
            C = connectionmatrix(femm1, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            np = Int(round(length(nl) / psize))
            # This partitions all the nodes, not just connected by the elements
            # above. 
            ppartitioning = Metis.partition(g, np; alg=:KWAY)
            # Therefore, we need to relabel the partition numbers. All the
            # partitions for nodes not connected together by fes will be set to
            # zero.
            epts = ppartitioning[nl]
            ppartitioning = fill(0, count(fens))
            ppartitioning[nl] = epts
            # Now we have the the current partitions. They are not necessarily
            # numbered consecutively. That needs to be fixed now.
            ppartitions = unique(epts)
            newppartitions = fill(0, maximum(ppartitions))
            newppartitions[ppartitions] = 1:length(ppartitions)
            for k in eachindex(ppartitioning)
                if ppartitioning[k] > 0
                    ppartitioning[k] = newppartitions[ppartitioning[k]]
                end 
            end 
            ppartitions = unique(ppartitioning[nl])
        end 
        # These are now the partitions of the nodes in the current domain.
        nap = length(allpartitions) + length(ppartitions)
        allpartitions = vcat(allpartitions, ppartitions)
        @assert nap == length(allpartitions)
        for k in eachindex(nl)
            gpartitioning[nl[k]] = ppartitioning[nl[k]] + cpartoffset
        end 
        cpartoffset += length(ppartitions)
    end 

    # Now we have the global partitioning of the nodes. Check that all nodes are
    # included (i.e. the partition is not 0)
    @assert isempty(findall(v -> v == 0, gpartitioning))
    # Some partitions may be missing: they are not necessarily in sequential
    # order. We need to remember them.
    gpartitions = unique(gpartitioning)
    newgpartitions = fill(0, maximum(gpartitions))
    newgpartitions[gpartitions] = 1:length(gpartitions)
    for k in eachindex(gpartitioning)
        if gpartitioning[k] > 0
            gpartitioning[k] = newgpartitions[gpartitioning[k]]
        end 
    end 
    gpartitions = unique(gpartitioning)
    # @show gpartitioning

    # that now we will start using the nodal partitioning.
    for k in eachindex(gpartitioning)
        @assert gpartitioning[k] in gpartitions
    end 

    X = fens.xyz
    mor = CoNCData(X, gpartitioning)

    P = NodalField(zeros(size(X, 1), 1)) # Pressure field
    numberdofs!(P) #
    Phi = transfmatrix(mor, LegendreBasis, n -> n == 1 ? (1:1) : (1:2), P)

    for j in 1:size(X, 1)
        P.values[j] = pressure(X[j, :])
    end 
    vtkexportmesh("rblock-Y.vtk", fens, fes; scalars = [("P", deepcopy(P.values))])

    Pa = Phi * ((Phi' * Phi) \ (Phi' * P.values))
    vtkexportmesh("rblock-Ya.vtk", fens, fes; scalars = [("Pa", deepcopy(Pa))])


    for j in 1:length(gpartitions)
        l = findall(v -> v == gpartitions[j], gpartitioning)
        sfes = FESetP1(reshape(l, length(l), 1))
        vtkexportmesh("rblock-np$(gpartitions[j]).vtk", fens, sfes)
    end


end
test()
end


module mesh_Q4multi_5
using Metis
using FinEtools
using FinEtools.MeshExportModule
using CoNCMOR: CoNCData, transfmatrix, LegendreBasis
using Statistics
using LinearAlgebra
using SparseArrays
using PlotlyLight
using Test

function pressure(xy)
    sin(7.5*xy[1]) + sin(5.3*xy[2])
end

function test()
    A, B = 2.0, 3.0
    nA, nB = 51, 73
    nelgroups = 5
    psize = 50
    dontpartition = [2, 4]

    fens, fes = Q4block(A, B, nA, nB)

    # Here we create the domain partitioning based on the finite elements.
    femm1 = FEMMBase(IntegDomain(fes, GaussRule(2, 1)))
    C = dualconnectionmatrix(femm1, fens, 2)
    g = Metis.graph(C; check_hermitian=true)
    dpartitioning = Metis.partition(g, nelgroups; alg=:KWAY)
    dpartitions = unique(dpartitioning)

    for j in eachindex(dpartitions)
        sfes = subset(fes, findall(v -> v == dpartitions[j], dpartitioning))
        vtkexportmesh("rblock-domain$(dpartitions[j]).vtk", fens, sfes)
    end

    gpartitioning = fill(0, count(fens))
    allpartitions = Int[]
    cpartoffset = 0
    for j in eachindex(dpartitions)
        sfes = subset(fes, findall(v -> v == dpartitions[j], dpartitioning))
        nl = connectednodes(sfes)
        if dpartitions[j] in dontpartition
            ppartitioning = fill(0, count(fens))
            ppartitions = 1:length(nl)
            ppartitioning[nl] .= ppartitions
        else
            femm1 = FEMMBase(IntegDomain(sfes, GaussRule(2, 1)))
            C = connectionmatrix(femm1, count(fens))
            g = Metis.graph(C; check_hermitian=true)
            np = Int(round(length(nl) / psize))
            # This partitions all the nodes, not just connected by the elements
            # above.
            ppartitioning = Metis.partition(g, np; alg=:KWAY)
            # Therefore, we need to relabel the partition numbers. All the
            # partitions for nodes not connected together by fes will be set to
            # zero.
            epts = ppartitioning[nl]
            ppartitioning = fill(0, count(fens))
            ppartitioning[nl] = epts
            # Now we have the the current partitions. They are not necessarily
            # numbered consecutively. That needs to be fixed now.
            ppartitions = unique(epts)
            newppartitions = fill(0, maximum(ppartitions))
            newppartitions[ppartitions] = 1:length(ppartitions)
            for k in eachindex(ppartitioning)
                if ppartitioning[k] > 0
                    ppartitioning[k] = newppartitions[ppartitioning[k]]
                end
            end
            ppartitions = unique(ppartitioning[nl])
        end
        # These are now the partitions of the nodes in the current domain.
        nap = length(allpartitions) + length(ppartitions)
        allpartitions = vcat(allpartitions, ppartitions)
        @assert nap == length(allpartitions)
        for k in eachindex(nl)
            gpartitioning[nl[k]] = ppartitioning[nl[k]] + cpartoffset
        end
        cpartoffset += length(ppartitions)
    end

    # Now we have the global partitioning of the nodes. Check that all nodes are
    # included (i.e. the partition is not 0)
    @assert isempty(findall(v -> v == 0, gpartitioning))
    # Some partitions may be missing: they are not necessarily in sequential
    # order. We need to renumber them.
    gpartitions = unique(gpartitioning)
    newgpartitions = fill(0, maximum(gpartitions))
    newgpartitions[gpartitions] = 1:length(gpartitions)
    for k in eachindex(gpartitioning)
        if gpartitioning[k] > 0
            gpartitioning[k] = newgpartitions[gpartitioning[k]]
        end
    end
    gpartitions = unique(gpartitioning)
    # @show gpartitioning

    # that now we will start using the nodal partitioning.
    for k in eachindex(gpartitioning)
        @assert gpartitioning[k] in gpartitions
    end

    X = fens.xyz
    mor = CoNCData(X, gpartitioning)

    P = NodalField(zeros(size(X, 1), 1)) # Pressure field
    numberdofs!(P) #
    Phi = transfmatrix(mor, LegendreBasis, n -> n == 1 ? (1:1) : (1:3), P)

    I, J, V = findnz(Phi)
    p = PlotlyLight.Plot()
    p(x = J, y = I, mode="markers")
    # p.layout.title.text = "Matrix G"
    p.layout.yaxis.title = "Row"
    p.layout.yaxis.range = [size(Phi, 1)+1, 0]
    p.layout.xaxis.title = "Column"
    p.layout.xaxis.range = [0, size(Phi, 2)+1]
    p.layout.xaxis.side = "top"
    p.layout.margin.pad = 10
    display(p)

    for j in 1:size(X, 1)
        P.values[j] = pressure(X[j, :])
    end
    vtkexportmesh("rblock-Y.vtk", fens, fes; scalars = [("P", deepcopy(P.values))])

    Pa = Phi * ((Phi' * Phi) \ (Phi' * P.values))
    vtkexportmesh("rblock-Ya.vtk", fens, fes; scalars = [("Pa", deepcopy(Pa))])


    # for j in 1:length(gpartitions)
    #     l = findall(v -> v == gpartitions[j], gpartitioning)
    #     sfes = FESetP1(reshape(l, length(l), 1))
    #     vtkexportmesh("rblock-np$(gpartitions[j]).vtk", fens, sfes)
    # end


end
test()
end

nothing
