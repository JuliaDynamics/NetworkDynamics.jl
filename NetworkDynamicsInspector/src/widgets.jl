struct TimeSlider{T} <: Bonito.AbstractSlider{T}
    values::Observable{Vector{T}}
    value::Observable{T}
    style::Styles
    track_style::Styles
    thumb_style::Styles
    track_active_style::Styles
end

function TimeSlider(value::Observable{T}, values::Observable{Vector{T}}) where {T}
    slider_height=15
    thumb_width=slider_height
    thumb_height=slider_height
    track_height=slider_height / 2
    track_active_height=track_height + 2
    backgroundcolor="transparent"
    track_color="#eee"
    track_active_color="#ddd"
    thumb_color="#fff"

    half_thumb_width = thumb_width / 2

    style = Styles(
        "display" => "grid",
        "grid-template-columns" => "1fr",
        "grid-template-rows" => "$(slider_height)px",
        "align-items" => "center",
        "margin" => "$(slider_height ÷ 3)px",
        "position" => "relative",
        "padding-right" => "$(2 + half_thumb_width)px",
        "padding-left" => "$(2 + half_thumb_width)px",
        "background-color" => backgroundcolor,
    )

    track_style = Styles(
        "position" => "absolute",
        "width" => "100%",
        "height" => "$(track_height)px",
        "background-color" => track_color,
        "border-radius" => "3px",
        "border" => "1px solid #ccc",
    )

    track_active_style = Styles(
        "position" => "absolute",
        "width" => "0px",
        "height" => "$(track_active_height)px",
        "background-color" => track_active_color,
        "border-radius" => "3px",
        "border" => "1px solid #ccc",
    )

    thumb_style = Styles(
        "width" => "$(thumb_width)px",
        "height" => "$(thumb_height)px",
        "background-color" => "white",
        "border-radius" => "50%",
        "cursor" => "pointer",
        "position" => "absolute",
        "border" => "1px solid #ccc",
        "left" => "$(-half_thumb_width)px",
        "background-color" => thumb_color,
    )

    slider = TimeSlider(
        values,
        value,
        style,
        track_style,
        thumb_style,
        track_active_style,
    )
    return slider
end

function Bonito.jsrender(session::Session, slider::TimeSlider)
    # Define the CSS styles
    container_style = slider.style
    track_style = slider.track_style
    track_active_style = slider.track_active_style
    thumb_style = slider.thumb_style

    # Create elements
    thumb = DOM.div(; style=thumb_style)
    track = DOM.div(; style=track_style)
    track_active = DOM.div(; style=track_active_style)
    container = DOM.div(track, track_active, thumb; style=container_style)

    # JavaScript for interactivity
    jscode = js"""
    (container)=> {
        const thumb = $(thumb);
        const track_active = $(track_active);
        const track = $(track);
        let isDragging = false;
        let thumbpos = 0;
        let thumbval = 0;
        $(slider.value).on(val => set_thumb_val(val));
        function move_thumb(e) {
            const width = track.offsetWidth;
            let new_left = e.clientX - container.getBoundingClientRect().left;
            new_left = Math.max(new_left, 0);
            new_left = Math.min(new_left, width);
            set_thumb_pos(new_left)
        }
        function set_thumb_pos(new_left) {
            if(new_left !== thumbpos) {
                const thumb_width = thumb.offsetWidth / 2;
                const width = track.offsetWidth;
                thumb.style.left = (new_left - thumb_width) + 'px';  // Update the left position of the thumb
                track_active.style.width = new_left + 'px';  // Update the active track

                const startval = $(slider.values[][begin]);
                const endval = $(slider.values[][end]);
                const new_val = startval + (new_left / width) * (endval - startval);

                thumbpos = new_left;
                thumbval = new_val;
                $(slider.value).notify(new_val);
            }
        }
        function set_thumb_val(new_val) {
            if(new_val !== thumbval) {
                const startval = $(slider.values[][begin]);
                const endval = $(slider.values[][end]);
                const thumb_width = thumb.offsetWidth / 2;
                const width = track.offsetWidth;
                const new_left = (new_val - startval) / (endval - startval) * width;
                thumb.style.left = (new_left - thumb_width) + 'px';  // Update the left position of the thumb
                track_active.style.width = new_left + 'px';  // Update the active track

                thumbpos = new_left;
                thumbval = new_val;
            }
        }
        function move_thumb_incremental(direction, shiftPressed) {
            const startval = $(slider.values[][begin]);
            const endval = $(slider.values[][end]);
            let step = (endval - startval) / 100;
            if (shiftPressed) {step *= 10;}
            if (direction === 'right') {
                set_thumb_val(Math.min(thumbval + step, endval));
            } else if (direction === 'left') {
                set_thumb_val(Math.max(thumbval - step, startval));
            }
        }
        const move_thumb_throttled = Bonito.throttle_function(move_thumb, 100);
        const controller = new AbortController();
        document.addEventListener('mousedown', function (e) {
            if(e.target === thumb || e.target === track_active || e.target === track || e.targar === container){
                isDragging = true;
                move_thumb(e);
                e.preventDefault();  // Prevent default behavior
            }
        }, { signal: controller.signal});
        document.addEventListener('mouseup', function () {
            if (!document.body.contains(container)) {
                controller.abort();
            }
            isDragging = false;
        }, { signal: controller.signal });
        document.addEventListener('mousemove', function (e) {
            if (isDragging) {
                move_thumb_throttled(e);
            }
        }, { signal: controller.signal });
        document.addEventListener('keydown', function (e) {
            if (e.key === 'ArrowRight') {
                move_thumb_incremental('right', e.shiftKey);
                e.preventDefault();
            } else if (e.key === 'ArrowLeft') {
                move_thumb_incremental('left', e.shiftKey);
                e.preventDefault();
            }
        }, { signal: controller.signal });
        set_thumb_val($(slider.value).value)
    }
    """
    onload(session, container, jscode)

    return jsrender(session, container)
end
