module ElectronExt

using NetworkDynamicsInspector: NetworkDynamicsInspector as NDI
using Electron: Electron, windows
using Bonito: Bonito, HTTPServer

const ELECTRON_APP = Ref{Any}(nothing)
const ELECTRON_DISP = Ref{Any}(nothing)

function NDI.serve_app(::NDI.ElectronDisp, app)
    disp = get_electron_display()
    display(disp, app)
    nothing
end

function NDI.close_display(::NDI.ElectronDisp; strict)
    if haswindow()
        @info "Close existing Windows"
        while length(windows(ELECTRON_APP[]))>0
            close(first(windows(ELECTRON_APP[])))
        end
    end
    if strict
        close_application()
    end
end

function get_electron_display()
    window = NDI.get_electron_window()
    # BUG: Electron display cannot be reused
    # if isnothing(ELECTRON_DISP[]) || window != ELECTRON_DISP[].window
    #     disp = HTTPServer.ElectronDisplay(window, HTTPServer.BrowserDisplay(; open_browser=false))
    #     ELECTRON_DISP[] = disp
    # else
    #     ELECTRON_DISP[]
    # end
    return HTTPServer.ElectronDisplay(window, HTTPServer.BrowserDisplay(; open_browser=false))
end

function NDI.get_electron_window()
    app = NDI.get_electron_app()

    any(w -> !w.exists, windows(app)) && @warn "App contains reference to nonexistent window(s)"

    window = if isempty(windows(app))
        x, y = NDI.CURRENT_DISPLAY[] isa NDI.ElectronDisp ? NDI.CURRENT_DISPLAY[].resolution : (1200, 800)
        opts = Dict(:width => x, :height => y, :webPreferences => Dict(:enableRemoteModule => true))
        @info "Create new Electron Window with $opts"
        Electron.Window(app, opts)
    else
        length(windows(app)) != 1 && @warn "App contains multiple windows"
        first(windows(app))
    end

    return window
end
haswindow() = hasapp() && !isempty(windows(ELECTRON_APP[]))

function NDI.get_electron_app()
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

function close_application()
    if hasapp()
        @info "Close Electron Application"
        close(ELECTRON_APP[])
    end
end

function NDI.toggle_devtools()
    if haswindow()
        Electron.toggle_devtools(NDI.get_electron_window())
    else
        error("No window to toggle devtools!")
    end
end
NDI.has_electron_window() = haswindow()


end
