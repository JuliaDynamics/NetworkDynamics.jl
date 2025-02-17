using Bonito
using NetworkDynamicsInspector
using NetworkDynamicsInspector: OptionGroup, MultiSelect

gui = (;
    range = Observable{NTuple{2, Float64}}((-1.0, 1.0)),
    val1 = Observable{Float64}(-1.0),
    val2 = Observable{Float64}(1.0),
)

let
    app = App(;) do session
        NetworkDynamicsInspector.clear_obs!(gui)
        slider = ContinuousSlider(gui.range, gui.val1, gui.val2)

        return wrap_assets(Grid(RoundedLabel(gui.val1), slider, RoundedLabel(gui.val2);
            columns="1fr 3fr 1fr", width="500px"))
    end;
    serve_app(app)
end

gui.val1[] = -0.5
gui.val2[] =  0.5
gui.val1[] = -1
gui.val2[] =  1

gui.range[] = (-0.5, 0.5)
gui.range[] = (-1, 1)
