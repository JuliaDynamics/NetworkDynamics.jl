using Bonito
using Electron
using NetworkDynamicsInspector
using NetworkDynamicsInspector: ContinuousSlider, RoundedLabel

function electron(app)
    close.(Electron.applications())
    disp = Bonito.use_electron_display(devtools = true)
    display(disp, app)
end

@testset "Test Slider" begin
    gui = (;
           range = Observable{NTuple{2, Float64}}((-1.0, 1.0)),
           val1 = Observable{Float64}(-1.0),
           val2 = Observable{Float64}(1.0),
           val = Observable{Float64}(0.0)
    )

    App(;) do session
        NetworkDynamicsInspector.clear_obs!(gui)
        sl1 = ContinuousSlider(gui.range, gui.val1, gui.val2)
        sl2 = ContinuousSlider(gui.range, gui.val; arrowkeys=true)

        Bonito.Grid(
            RoundedLabel(gui.val1), sl1, RoundedLabel(gui.val2),
            RoundedLabel(gui.val), sl2;
            columns="1fr 3fr 1fr", width="500px"
        )
    end |> electron

    gui.val1[] = -0.5
    gui.val2[] =  0.5
    gui.val1[] = -1
    gui.val2[] =  1

    gui.range[] = (-0.5, 0.5)
    gui.range[] = (-1, 1)
end
