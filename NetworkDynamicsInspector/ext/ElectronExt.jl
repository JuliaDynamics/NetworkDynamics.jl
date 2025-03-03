module ElectronExt

using NetworkDynamicsInspector: NetworkDynamicsInspector as NDI
using Electron: Electron
using Bonito: Bonito

function NDI.display_electron_app()
    webapp = NDI.get_webapp()
    # disp = Bonito.use_electron_display()
    disp = Bonito.use_electron_display(devtools = true)
    display(disp, webapp)
    nothing
end

end
