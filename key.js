
sketch.default2d();

// Define colors
var white_off_color = [1., 1., 1.];
var white_on_color = [0.4, 0.6, 0.4];
var black_off_color = [0., 0., 0.];
var black_on_color = [0.2, 0.4, 0.2];

// Initialize zeros and click coordinates
var octave = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.];
var is_black_key = [false, true, false, true, false, false, true, false, true, false, true, false];
var last_x_world = 0;
var last_y_world = 0;

var num_keys = 12;

var width = box.rect[2] - box.rect[0];
var height = box.rect[3] - box.rect[1];

draw();

// Toggle the value on double click
function ondblclick(x_screen, y_screen) {
    var key = get_key(x_screen, y_screen);
    if (octave[key] == 1.) {
        octave[key] = 0.;
    } else {
        octave[key] = 1.;
    }
    onclick(x_screen, y_screen);
    bang();
}

// Store the last click to track deltas
function onclick(x_screen, y_screen) {
    last_key = get_key(x_screen, y_screen);
	last_x_screen = x_screen;
	last_y_screen = y_screen;
}

// Change the clicked value based on drag
function ondrag(x_screen, y_screen) {
	// calculate delta movements
	var dy = y_screen - last_y_screen;
    
    octave[last_key] = octave[last_key] - dy * 0.005;
    octave[last_key] = Math.min(1, Math.max(0, octave[last_key]))

	last_x_screen = x_screen;
	last_y_screen = y_screen;
    bang();
}

function screen_to_relative(x_screen, y_screen) {
    return [x_screen/width, y_screen/height];
}

function relative_to_world(x, y) {
    var x_relative = (2 * x - 1) * width/height;
    var y_relative = 2 * y - 1;
    return [x_relative, y_relative];
}

// Get the key at the specified coordinates
function get_key(x_screen, y_screen) {
    var x = screen_to_relative(x_screen, y_screen)[0];
    var key = Math.floor(x * num_keys);
    key = Math.min(num_keys, Math.max(0, key));
    return key;
}

function draw() {
	with (sketch) {

        // Clear the screen
		glclear();

        for (var key = 0; key < num_keys; key++) {
            var color = [];
            if (is_black_key[key]) {
                var on_color = black_on_color;
                var off_color = black_off_color;
            } else {
                var on_color = white_on_color;
                var off_color = white_off_color;
            }
            for(var i = 0; i < 3; i++) {
                color.push(on_color[i] * octave[key] + off_color[i] * (1 - octave[key]));
            }
            glcolor(color);

            xy_min = relative_to_world(key/num_keys, 0);
            xy_max = relative_to_world((key + 1)/num_keys, 1);

            
            glrect(xy_min[0], xy_min[1], xy_max[0], xy_max[1]);
        }
	}
}

function list(x) {
    if (arguments.length != num_keys) return;
    for (var key = 0; key < num_keys; key++) {
        octave[key] = arguments[key];
    }
    bang();
}

function bang() {
	draw();
	refresh();
	outlet(0,octave);
}

function onresize(w,h) {
    width = box.rect[2] - box.rect[0];
	height = box.rect[3] - box.rect[1];
	draw();
	refresh();
}