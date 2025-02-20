using NetworkDynamicsInspector, Bonito

tog = Observable{Bool}(true)
let
    _app = App() do session
        NetworkDynamicsInspector.clear_obs!(tog)
        toggle = NetworkDynamicsInspector.ToggleSwitch(value=tog, label="Toggle me")
        on(toggle.value) do state
            @info "value = $state"
        end
        DOM.div(
            DOM.span("text before"),
            toggle,
        )
    end;
    serve_app(_app)
end
tog[] = !tog[]
