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
    arrowkeys::Bool
end

function ContinuousSlider(range, value::Observable{T}; kwargs...) where {T}
    value_l = Observable{T}(zero(T)/zero(T))
    ContinuousSlider(range, value_l, value; kwargs...)
end

function ContinuousSlider(range, value_l::Observable{T}, value_r::Observable{T}; arrowkeys=false) where {T}
    style, track_style, track_active_style, thumb_style = get_slider_style()

    #=
    on(range; update=true) do _r
        changed = false
        l_c = clamp(value_l[], _r[1], _r[2])
        if l_c != value_l
            changed = true
            value_l[] = l_c
        end
        r_c = clamp(value_r[], _r[1], _r[2])
        if r_c != value_r
            changed = true
            value_r[] = r_c
        end
        changed && @debug "Slider: range changed => clamped thumb values"
        nothing
    end
    =#

    slider = ContinuousSlider(
        range,
        value_l,
        value_r,
        style,
        track_style,
        thumb_style,
        track_active_style,
        arrowkeys
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

        function move_thumb_incremental(direction, shiftPressed) {
            const startval = $(slider.range).value[0];
            const endval = $(slider.range).value[1];
            let step = (endval - startval) / 100;
            if (shiftPressed) {step *= 10;}

            let newval;
            if (direction === 'right') {
                newval = Math.min(thumbval_r + step, endval);
            } else if (direction === 'left') {
                newval = Math.max(thumbval_r - step, startval);
            }
            $(slider.value_r).notify(newval);
        }

        if ($(slider.arrowkeys)) {
            document.addEventListener('keydown', function (e) {
                if (e.key === 'ArrowRight') {
                    move_thumb_incremental('right', e.shiftKey);
                    e.preventDefault();
                } else if (e.key === 'ArrowLeft') {
                    move_thumb_incremental('left', e.shiftKey);
                    e.preventDefault();
                }
            }, { signal: controller.signal });
        }

        // update positions on resize
        const observer = new ResizeObserver(entries => {
            thumbpos_r = isNaN(thumbval_r) ? 0 : val_to_pos(thumbval_r);
            thumbpos_l = isNaN(thumbval_l) ? 0 : val_to_pos(thumbval_l);

            const thumb_width = thumb_r.offsetWidth / 2;
            const left_r = (thumbpos_r - thumb_width) + 'px';
            const left_l = (thumbpos_l - thumb_width) + 'px';

            thumb_l.style.left = left_l;
            thumb_r.style.left = left_r;
            track_active.style.left = thumbpos_l + 'px';  // Update the active track
            track_active.style.width = (thumbpos_r-thumbpos_l) + 'px';  // Update the active track
        });
        observer.observe(container);
    }
    """
    onload(session, container, jscode)

    return jsrender(session, container)
end


function RoundedLabel(value; sigdigits=2, style=Styles(), attributes...)
    styled = Styles(style, "font-size" => "1rem")
    str = @lift NetworkDynamics.str_significant($value; sigdigits=$sigdigits)
    return DOM.span(str; style=styled)
end

struct OptionGroup{T}
    label::String
    options::Vector{T}
end

struct MultiSelect{T}
    options::Observable{Vector{Union{T, OptionGroup{T}}}}
    selection::Observable{Vector{T}}
    placeholder::String
    multi::Bool
    option_to_string::Any
    function MultiSelect(_options, _selection=nothing; T=Any, placeholder="", multi=true, option_to_string=repr)
        options = _options isa Observable ? _options : Observable{Vector{T}}(_options)
        selection = if isnothing(_selection)
            Observable(T[])
        elseif _selection isa Observable
            @assert _selection isa Observable{Vector{T}}
            _selection
        else
            Observable{Vector{T}}(_selection)
        end
        new{T}(options, selection, placeholder, multi, option_to_string)
    end
end

function Bonito.jsrender(session::Session, multiselect::MultiSelect)
    # jquery = Asset("https://cdn.jsdelivr.net/npm/jquery@3.7.1/dist/jquery.min.js")
    # select2_css = Asset("https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/css/select2.min.css")
    # select2_js = Asset("https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/js/select2.min.js")

    # generate internal observables of js representation of options and selection
    jsoptions = @lift options_to_jsoptions($(multiselect.options); option_to_string=multiselect.option_to_string)
    jsselection = Observable{Vector{Int}}()

    onany(jsselection) do _jssel
        sel = jsselection_to_selection(multiselect.options[], _jssel)
        # @info "New jsselection triggers new selection:" _jssel sel
        if sel != multiselect.selection[]
            multiselect.selection[] = sel
        end
        nothing
    end

    onany(multiselect.options, multiselect.selection; update=true) do _opts, _sel
        if !multiselect.multi && length(_sel) > 1
            deleteat!(_sel, firstindex(_sel)+1:lastindex(_sel))
        end
        jssel = selection_to_jsselection(_opts, _sel)
        @debug "MS \"$(multiselect.placeholder)\": New opts or sel triggers new jsselection:" _opts _sel jssel

        # filter out invalid selections
        if any(isnothing, jssel)
            invalid = findall(isnothing, jssel)
            deleteat!(_sel, invalid)
            filter!(!isnothing, jssel)
        end

        jsselection[] = jssel
        nothing
    end

    # Create a multi-select element
    id = replace(string(gensym("selectbox")), "#"=>"")

    style = Styles(
        "width" => "100%",
        # "font-family" => "monospace",
    )

    class = multiselect.multi ? "bonito-select2-multi" : "bonito-select2-single"
    select = DOM.select(
        DOM.option();
        multiple = multiselect.multi,
        class,
        # name = "states[]",
        style,
        id
    )

    container = DOM.div(
        # jquery,
        # select2_css,
        # select2_js,
        select
    )

    jqdocument = Bonito.JSString(raw"$(document)")
    jqselect = Bonito.JSString(raw"$('#"* id * raw"')")

    esc = Bonito.JSString(raw"$")
    js_onload = js"""
    (select) => {
        $(jqdocument).ready(function(){
            $jqselect.select2({
                width: '100%',
                placeholder: $(multiselect.placeholder)
            });
        });

        function array_equal(a1, a2){
            return a1.length === a2.length &&
                    a1.every(function(val, i) { return val == a2[i]})
        }

        function updateDisplayedOptions(new_options) {
            const jq_select = $jqselect;

            // Clear previous options
            jq_select.empty();

            // Loop through each group and create optgroups
            new_options.forEach(opt_or_group => {
                if (opt_or_group.options === undefined) {
                    let newOption = new Option(opt_or_group.text, opt_or_group.id, false, false);
                    jq_select.append(newOption);
                } else {
                    let jq_optgroup = $esc('<optgroup>', { label: opt_or_group.label });

                    opt_or_group.options.forEach(option => {
                        let newOption = new Option(option.text, option.id, false, false);
                        jq_optgroup.append(newOption);
                    });

                    jq_select.append(jq_optgroup);
                }
            });
        }
        updateDisplayedOptions($(jsoptions).value);
        $(jsoptions).on(updateDisplayedOptions);

        function updateDisplayedSelection(new_sel_nr) {
            const jq_select = $jqselect
            const new_sel = new_sel_nr.map(String);

            jq_select.data('preserved-order', new_sel);
            jq_select.val(new_sel).trigger('change');
            select2_renderSelections();
        }
        updateDisplayedSelection($(jsselection).value)
        $(jsselection).on(updateDisplayedSelection)

        // Trigger update for Select2 to recognize new options
        $jqselect.trigger('change');

        // Don't reorder
        // https://github.com/select2/select2/issues/3106#issuecomment-333341636
        function select2_renderSelections(){
            jq_select2 = $jqselect
            const def_order  = jq_select2.val();
            const pre_order  = jq_select2.data('preserved-order') || [];
            const jq_tags    = jq_select2.next('.select2-container').find('li.select2-selection__choice');
            const jq_tags_ul = jq_tags.first().parent()

            const new_order = pre_order.map(val=>{
                return def_order.indexOf(val);
            });

            const sortedElements = new_order.map(i => jq_tags.eq(i));
            jq_tags_ul.append(sortedElements);
        }
        function selectionHandler(e){
            const jq_select2  = $esc(this);
            const val         = e.params.data.id;

            let order;
            if ($(multiselect.multi)){
                order = jq_select2.data('preserved-order') || [];
            } else {
                order = [];
            }

            switch (e.type){
                case 'select2:select':
                    order[ order.length ] = val;
                    break;
                case 'select2:unselect':
                    let found_index = order.indexOf(val);
                    if (found_index >= 0 )
                        order.splice(found_index,1);
                    break;
            }

            jq_select2.data('preserved-order', order); // store it for later
            select2_renderSelections();

            // notify julia about changed selection
            $(jsselection).notify(order.map(Number));
        }
        $jqselect.on('select2:select select2:unselect', selectionHandler);
    }
    """
    Bonito.onload(session, select, js_onload)

    return jsrender(session, container)
end

function options_to_jsoptions(options; option_to_string=repr)
    jsoptions = []
    id = 1
    for option in options
        if option isa OptionGroup
            jssuboptions = []
            for suboption in option.options
                push!(jssuboptions, (;text=option_to_string(suboption), id=id))
                id += 1
            end
            push!(jsoptions, (;label=option.label, options=jssuboptions))
        else
            push!(jsoptions, (;text=option_to_string(option), id=id))
            id += 1
        end
    end
    jsoptions
end

jsselection_to_selection(options, jsselection) = _jsselection_to_selection.(Ref(options), jsselection)
function _jsselection_to_selection(options, jsselection::Int)
    id = 1
    for option in options
        if option isa OptionGroup
            for suboption in option.options
                id == jsselection && return suboption
                id += 1
            end
        else
            id == jsselection && return option
            id += 1
        end
    end
end

selection_to_jsselection(options, selection) = _selection_to_jsselection.(Ref(options), selection)
function _selection_to_jsselection(options, selection)
    id = 1
    for option in options
        if option isa OptionGroup
            for suboption in option.options
                suboption == selection && return id
                id += 1
            end
        else
            option == selection && return id
            id += 1
        end
    end
end
