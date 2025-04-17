using NetworkDynamicsInspector, Bonito

tog = Observable{Bool}(true)
let
    _app = App() do session
        NetworkDynamicsInspector.clear_obs_and_close!(tog)
        toggle = NetworkDynamicsInspector.ToggleSwitch(value=tog, label="Toggle me")
        on(toggle.value) do state
            @info "value = $state"
        end
        toggle2 = NetworkDynamicsInspector.ToggleSwitch(label="long labeel")
        DOM.div(
            DOM.span("text before"),
            toggle,
            DOM.div(toggle2, style=Styles("width"=>"100px"))
        )
    end;
    serve_app(_app)
end
tog[] = !tog[]
