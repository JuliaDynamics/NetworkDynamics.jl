function isElectron() {
    // Check for the presence of Electron-specific objects
    return typeof window !== 'undefined' && window.process && window.process.type === 'renderer';
}

// Function to handle zoom in/out
function handleZoom(event) {
    if (event.ctrlKey) {
        if (event.key === '+' || event.key === '=') {
            // Zoom in
            const { webFrame } = require('electron');
            webFrame.setZoomLevel(webFrame.getZoomLevel() + 1);
            event.preventDefault();
        } else if (event.key === '-' || event.key === '_') {
            // Zoom out
            const { webFrame } = require('electron');
            webFrame.setZoomLevel(webFrame.getZoomLevel() - 1);
            event.preventDefault();
        }
    }
}

// Function to open the developer console
function openDevTools(event) {
    if (event.ctrlKey && event.shiftKey && event.key === 'K') {
        console.log("Opening developer tools");
        const { getCurrentWindow } = require('@electron/remote');
        const currentWindow = getCurrentWindow();
        currentWindow.webContents.toggleDevTools();
        event.preventDefault();
    }
}

// Attach event listeners only if running in Electron
if (isElectron()) {
    document.addEventListener('keydown', handleZoom);
    //document.addEventListener('keydown', openDevTools);
}
