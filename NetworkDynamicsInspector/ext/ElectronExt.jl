module ElectronExt

using NetworkDynamicsInspector: NetworkDynamicsInspector as NDI
using Electron: Electron, windows
using Bonito: Bonito, HTTPServer

const ELECTRON_APP = Ref{Any}(nothing)

function NDI.display_electron_app(restart)
    restart && close_windows()
    webapp = NDI.get_webapp()
    disp = get_electron_display()
    display(disp, webapp)
    nothing
end

function get_electron_display()
    window = get_electron_window()
    disp = HTTPServer.ElectronDisplay(window, HTTPServer.BrowserDisplay(; open_browser=false))
end

function get_electron_window()
    app = get_electron_app()

    any(w -> !w.exists, windows(app)) && @warn "App contains reference to nonexistent window(s)"

    window = if isempty(windows(app))
        opts = Dict(:width => 1200, :height => 800)
        Electron.Window(app, opts)
    else
        length(windows(app)) != 1 && @warn "App contains multiple windows"
        first(windows(app))
    end
    return window
end
haswindow() = hasapp() && !isempty(windows(ELECTRON_APP[]))

function get_electron_app()
    if !hasapp()
        ELECTRON_APP[] = Electron.Application(;
            additional_electron_args=[
                "--disable-logging",
                "--no-sandbox",
                "--user-data-dir=$(mktempdir())",
                "--disable-features=AccessibilityObjectModel",
            ],
        )
    end
    ELECTRON_APP[]
end
hasapp() = !isnothing(ELECTRON_APP[]) && ELECTRON_APP[].exists

function close_windows()
    if haswindow()
        @info "Close existing Windows"
        close.(windows(ELECTRON_APP[]))
    end
end

function close_application()
    if hasapp()
        @info "Close Electron Application"
        close(ELECTRON_APP[])
    end
end

function NDI.toggle_devtools()
    if haswindow()
        Electron.toggle_devtools(get_electron_window())
    else
        error("No window to toggle devtools!")
    end
end

end
