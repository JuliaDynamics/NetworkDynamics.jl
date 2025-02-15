function get_slider_style()
    slider_height=15
    thumb_width=slider_height
    thumb_height=slider_height
    track_height=slider_height / 3
    track_active_height=track_height# + 2
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
    (; style, track_style, track_active_style, thumb_style)
end

struct ContinuousSlider{T} <: Bonito.AbstractSlider{T}
    range::Observable{Tuple{T,T}}
    value_l::Observable{T}
    value_r::Observable{T}
    style::Styles
    track_style::Styles
    thumb_style::Styles
    track_active_style::Styles
end

function ContinuousSlider(range, value::Observable{T}) where {T}
    value_l = Observable{T}(zero(T)/zero(T))
    ContinuousSlider(range, value_l, value)
end

function ContinuousSlider(range, value_l::Observable{T}, value_r::Observable{T}) where {T}
    style, track_style, track_active_style, thumb_style = get_slider_style()

    slider = ContinuousSlider(
        range,
        value_l,
        value_r,
        style,
        track_style,
        thumb_style,
        track_active_style,
    )
    return slider
end

function Bonito.jsrender(session::Session, slider::ContinuousSlider)
    # Define the CSS styles
    container_style = slider.style
    track_style = slider.track_style
    track_active_style = slider.track_active_style
    thumb_style = slider.thumb_style

    # Create elements
    thumb_l = DOM.div(; style=thumb_style)
    thumb_r = DOM.div(; style=thumb_style)
    track = DOM.div(; style=track_style)
    track_active = DOM.div(; style=track_active_style)
    container = DOM.div(track, track_active, thumb_l, thumb_r; style=container_style)

    # JavaScript for interactivity
    jscode = js"""
    (container)=> {
        const thumb_r = $(thumb_r);
        const thumb_l = $(thumb_l);
        const track_active = $(track_active);
        const track = $(track);
        let isDragging_r = false;
        let isDragging_l = false;
        let thumbpos_l = 0;
        let thumbval_l = 0;
        let thumbpos_r = 0;
        let thumbval_r = 0;
        $(slider.value_l).on(val => set_thumb_val(val, 'l'));
        $(slider.value_r).on(val => set_thumb_val(val, 'r'));
        $(slider.range).on(val => set_thumb_val(thumbval_l, 'l'));
        $(slider.range).on(val => set_thumb_val(thumbval_r, 'r'));

        function thumb_event(e, thumb) {
            let new_pos = e.clientX - container.getBoundingClientRect().left;
            const width = track.offsetWidth;
            if (thumb==='r') {
                llim = (e.shiftKey) ? (Math.ceil(width/2)+1) : thumbpos_l+1;
                rlim = width;
            } else if (thumb==='l') {
                llim = 0;
                rlim = (e.shiftKey) ? (Math.floor(width/2)-1) : thumbpos_r-1;
            }
            new_pos = Math.max(new_pos, llim);
            new_pos = Math.min(new_pos, rlim);
            set_thumb_pos(new_pos, thumb);

            if (e.shiftKey) {
                const other = thumb==='l' ? 'r' : 'l';
                set_thumb_pos(width - new_pos, other);
            }
        }

        function set_thumb_pos(new_pos, thumb) {
            if (thumb=='l' && new_pos !== thumbpos_l || thumb=='r' && new_pos !== thumbpos_r) {
                const thumb_width = thumb_r.offsetWidth / 2;
                const width = track.offsetWidth;
                const new_left = (new_pos - thumb_width) + 'px';

                const new_val = pos_to_val(new_pos)
                if (thumb=='l') {
                    thumbpos_l = new_pos;
                    thumbval_l = new_val;
                    thumb_l.style.left = new_left
                    $(slider.value_l).notify(new_val);
                } else {
                    thumbpos_r = new_pos;
                    thumbval_r = new_val;
                    thumb_r.style.left = new_left
                    $(slider.value_r).notify(new_val);
                }
                track_active.style.left = thumbpos_l + 'px';  // Update the active track
                track_active.style.width = (thumbpos_r-thumbpos_l) + 'px';  // Update the active track
            }
        }

        function set_thumb_val(new_val, thumb) {
            const thumb_width = thumb_r.offsetWidth / 2;
            const new_pos = isNaN(new_val) ? 0 : val_to_pos(new_val);
            const new_left = (new_pos - thumb_width) + 'px';
            if (thumb=='l') {
                if (isNaN(new_val)) {
                    thumb_l.style.display = 'none';
                } else if (isNaN(thumbval_l)) {
                    thumb_l.style.display = 'block';
                }
                thumbpos_l = new_pos;
                thumbval_l = new_val;
                thumb_l.style.left = new_left;
            } else {
                thumbpos_r = new_pos;
                thumbval_r = new_val;
                thumb_r.style.left = new_left;
            }
            track_active.style.left = thumbpos_l + 'px';  // Update the active track
            track_active.style.width = (thumbpos_r-thumbpos_l) + 'px';  // Update the active track
        }

        function pos_to_val(pos) {
            const width = track.offsetWidth;
            const startval = $(slider.range).value[0];
            const endval = $(slider.range).value[1];
            const val = startval + (pos / width) * (endval - startval);
            return val;
        }

        function val_to_pos(val) {
            const width = track.offsetWidth;
            const startval = $(slider.range).value[0];
            const endval = $(slider.range).value[1];
            const pos = (val - startval) / (endval - startval) * width;
            return pos;
        }

        //const move_thumb_throttled = Bonito.throttle_function(thumb_event, 100);
        const controller = new AbortController();
        document.addEventListener('mousedown', function (e) {
            if (e.target === thumb_r){
                isDragging_r = true;
                thumb_event(e, 'r');
                e.preventDefault();  // Prevent default behavior
            } else if (e.target === thumb_l) {
                isDragging_l = true;
                thumb_event(e, 'l');
                e.preventDefault();  // Prevent default behavior
            } else if (e.target == track_active || e.target === track){
                let new_pos = e.clientX - container.getBoundingClientRect().left;
                // if closer to the right thumb, move the right thumb
                if(isNaN(thumbval_l) || Math.abs(new_pos - thumbpos_r) < Math.abs(new_pos - thumbpos_l)) {
                    thumb_event(e, 'r');
                } else {
                    thumb_event(e, 'l');
                }
                e.preventDefault();  // Prevent default behavior
            }
        }, { signal: controller.signal});
        document.addEventListener('mouseup', function () {
            if (!document.body.contains(container)) {
                controller.abort();
            }
            isDragging_r = false;
            isDragging_l = false;
        }, { signal: controller.signal });
        document.addEventListener('mousemove', function (e) {
            if (isDragging_r) {
                thumb_event(e, 'r');
            } else if (isDragging_l) {
                thumb_event(e, 'l');
            }
        }, { signal: controller.signal });
        set_thumb_val($(slider.value_l).value, 'l');
        set_thumb_val($(slider.value_r).value, 'r');
    }
    """
    # function move_thumb_incremental(direction, shiftPressed) {
    #     const startval = $(slider.values[][begin]);
    #     const endval = $(slider.values[][end]);
    #     let step = (endval - startval) / 100;
    #     if (shiftPressed) {step *= 10;}
    #     if (direction === 'right') {
    #         set_thumb_val(Math.min(thumbval + step, endval));
    #     } else if (direction === 'left') {
    #         set_thumb_val(Math.max(thumbval - step, startval));
    #     }
    # }
    # document.addEventListener('keydown', function (e) {
    #     if (e.key === 'ArrowRight') {
    #         move_thumb_incremental('right', e.shiftKey);
    #         e.preventDefault();
    #     } else if (e.key === 'ArrowLeft') {
    #         move_thumb_incremental('left', e.shiftKey);
    #         e.preventDefault();
    #     }
    # }, { signal: controller.signal });
    onload(session, container, jscode)

    return jsrender(session, container)
end


function RoundedLabel(value; sigdigits=2, style=Styles(), attributes...)
    styled = Styles(style, "font-size" => "1rem")
    str = @lift NetworkDynamics.str_significant($value; sigdigits=$sigdigits)
    return DOM.span(str; style=styled)
end
