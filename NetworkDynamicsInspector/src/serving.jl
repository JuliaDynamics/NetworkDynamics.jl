abstract type NDIDisplay end

@kwdef struct ElectronDisp <: NDIDisplay
    resolution::Tuple{Int, Int} = (1200, 800)
end

struct BrowserDisp <: NDIDisplay end

@kwdef struct ServerDisp <: NDIDisplay
    port::Int = 8080
    url::String = "localhost"
end

function disp_from_sym(sym::Symbol)
    if sym == :Electron
        return ElectronDisp()
    elseif sym == :Browser
        return BrowserDisp()
    elseif sym == :Server
        return ServerDisp()
    else
        error("Unknown display type: $sym")
    end
end

"""
    get_display(::NDIDisplay; restart, kwargs...)
"""
function serve_app(display, app)
    if disp isa ElectronDisp
        @error "Electron.jl not available. Please install Electron.jl and `using Electron` before calling this function."
    else
        @error "Unknown display type: $disp"
    end
end
close_display(::Nothing; kwargs...) = nothing
close_display() = close_display(CURRENT_DISPLAY[]; strict=true)

####
#### Server Display
####
SERVER_STATE = Dict()
function serve_app(s::ServerDisp, app)
    server = get(SERVER_STATE, :server, nothing)
    served_app = get(SERVER_STATE, :app, nothing)

    if served_app !== app
        close_display(s)
    end

    new_server = false
    if isnothing(server) || !Bonito.HTTPServer.isrunning(server)
        new_server = true
        server = Bonito.Server(app, s.url, s.port)
    end

    SERVER_STATE[:server] = server
    SERVER_STATE[:app] = app
    if new_server
        @info "Inspector server started at http://$(s.url):$(s.port)"
    else
        @info "Inspector server already running at http://$(s.url):$(s.port)"
    end
end
function close_display(::ServerDisp; kwargs...)
    server = get(SERVER_STATE, :server, nothing)
    if !isnothing(server) && Bonito.HTTPServer.isrunning(server)
        close_session(SESSION[])
        @info "Close running server"
        close(server)
    end
end

####
#### Browser Display
####
BROWSER_STATE = Dict()
function serve_app(::BrowserDisp, app)
    handler = get(BROWSER_STATE, :handler, nothing)
    if isnothing(handler)
        handler = Bonito.HTTPServer.BrowserDisplay()
    end
    BROWSER_STATE[:handler] = handler
    display(handler, app)
    nothing
end
function close_display(::BrowserDisp; kwargs...)
    BROWSER_STATE[:handler] = nothing
end



function close_session(session)
    if !isnothing(session) && Base.isopen(session)
        @info "Close existing session"
        on_session_close = js"""
            // Create a semi-transparent overlay
            const overlay = document.createElement('div');
            overlay.style.position = 'fixed';
            overlay.style.top = '0';
            overlay.style.left = '0';
            overlay.style.width = '100%';
            overlay.style.height = '100%';
            overlay.style.backgroundColor = 'rgba(0, 0, 0, 0.7)'; // Semi-transparent black
            overlay.style.zIndex = '1000'; // Ensure it's on top of everything
            document.body.appendChild(overlay);

            // Create the message element
            const message = document.createElement('div');
            message.textContent = 'Session Closed';
            message.style.position = 'fixed';
            message.style.top = '50%';
            message.style.left = '50%';
            message.style.transform = 'translate(-50%, -50%)';
            message.style.color = 'white';
            message.style.fontSize = '2em';
            message.style.zIndex = '1001'; // Ensure it's on top of the overlay
            message.style.textAlign = 'center';
            overlay.appendChild(message);

            // Optionally, prevent any interaction with the underlying elements
            overlay.style.pointerEvents = 'auto';
        """
        Bonito.evaljs(session, on_session_close)
        close(session)
    end
end

# functions defined in extension
function toggle_devtools end
function has_electron_window end
function get_electron_app end
function get_electron_window end

function save_electron_screenshot(path=joinpath(@__DIR__, "screenshot.png"), resize=true)
    sync()
    if isempty(run(get_electron_window(),
        "Array.from(document.querySelectorAll('.graphplot-col, .timeseries-col'))"))
        error("No content to screenshot!")
    end
    resize && _resize_electron_to_content()

    path = isabspath(path) ? path : joinpath(pwd(), path)
    has_electron_window()|| error("No Electron window exists!")
    winid = get_electron_window().id
    js = """
    let win = BrowserWindow.fromId($winid)
    win.webContents.capturePage().then(image => {
        const screenshotPath = '$path';
        require('fs').writeFileSync(screenshotPath, image.toPNG());
        console.log('Screenshot saved to ', screenshotPath);
    });
    """
    sync()
    sleep(3) # make sure that axis updates and so on
    d = run(get_electron_app(), js)
    while !(isfile(path))
        sleep(0.1)
    end
    nothing
end

function _resize_electron_to_content()
    window = get_electron_window()
    js_max_h = """
    {
        let maxHeight = 0;
        Array.from(document.querySelectorAll('.graphplot-col, .timeseries-col')).forEach(el => {
            const height = el.offsetHeight;
            // console.log("Set new height", height)
            if (height > maxHeight) {
                maxHeight = height;
            }
        })
        maxHeight
    }
    """
    y = run(window, js_max_h)

    # set size
    app = get_electron_app()
    resize_js = """
    {
        let window = BrowserWindow.fromId($(window.id));
        let size = window.getSize();
        window.setSize(size[0], $y);
    }
    """
    run(app, resize_js)
    sleep(3)
end
