// Helper
function AMax(arr) {
  let result = arr[0];
  for (let i = 0; i < arr.length; i++)
      if (arr[i] > result)
          result = arr[i];
  return result;
}

function AMin(arr) {
  let result = arr[0];
  for (let i = 1; i < arr.length; i++)
      if (arr[i] < result)
          result = arr[i];
  return result;
}

function AMaxI(arr) {
  let result = arr[0];
  let imax = 0;
  for (let i = 1; i < arr.length; i++)
      if (arr[i] > result) {
          result = arr[i];
          imax = i;
        }
  return [result, imax];
}

function clr(arr, val) {
    for (let i = 0; i < arr.length; i++)
        arr[i] = val;
}

// The function gets called when the window is fully loaded
// Get the canvas and context
var canvas1 = document.getElementById("view1");
var context1 = canvas1.getContext("2d");
var canvas2 = document.getElementById("view2");
var context2 = canvas2.getContext("2d");

// Define the image dimensions
var width = canvas1.width;
var height = canvas1.height;

// Create an ImageData object
var imagedata1 = context1.createImageData(width, height);
var imagedata2 = context2.createImageData(width, height);

var n_cluster = 20;

var n_total = width * height;
var bias = 0;//n_cluster / n_total;//1/(width*height);

var pfield = new Uint32Array(width * height);
var prev_pfield = new Uint32Array(width * height);

var cfield = new Float64Array(width * height);
var prev_cfield = new Float64Array(width * height);

var pressure = new Array(n_cluster + 1);
var prev_pressure = new Array(n_cluster + 1);

var desired_count = new Array(n_cluster + 1);
var urf_c = 1;
var urf_p = 1;

function sin_col(i, size, phase) {
  var sin = Math.sin(Math.PI / (size / 2) * i + phase);
  var int = Math.floor(sin * 127) + 128;
  return int;
}

cb_p = new Array(n_cluster + 1);

var cb_c = new Array(256);
for (var i = 0; i < 256; i++) {
  var red   = Math.min(i/100,1) * sin_col(i, 255, 0 * Math.PI * 2/3); // 0   deg
  var blue  = Math.min(i/100,1) * sin_col(i, 255, 1 * Math.PI * 2/3); // 120 deg
  var green = Math.min(i/100,1) * sin_col(i, 255, 2 * Math.PI * 2/3); // 240 deg
  cb_c[i] = [red, green, blue];
}


// Init

// Percentage of domain
/*
pressure[1] = 0.5;
pressure[2] = 0.25;
pressure[3] = 0.125;
pressure[4] = 0.125/2;
pressure[5] = 0.125/2;
*/
function restart() {
  ts=0;
  n_cluster = Number(document.getElementById("np").value);
  desired_count = new Array(n_cluster + 1);
  cb_p = new Array(n_cluster + 1);
  for (var i = 0; i < n_cluster + 1; i++) {
    var red   = sin_col(i, n_cluster + 1, 0 * Math.PI * 2/3); // 0   deg
    var blue  = sin_col(i, n_cluster + 1, 1 * Math.PI * 2/3); // 120 deg
    var green = sin_col(i, n_cluster + 1, 2 * Math.PI * 2/3); // 240 deg
    cb_p[i] = [red, green, blue];
  }


  w0 = Number(document.getElementById("w0").value)
  w1 = Number(document.getElementById("w1").value)
  pressure = new Array(n_cluster + 1);
  prev_pressure = new Array(n_cluster + 1);
  weights = new Array(n_cluster);
  weights.fill(0);
  weights[0] = w0;
  weights[1] = w1;
  // auto balance
  let weights_sum = weights.reduce((ts, el)=>ts+=el, 0);
  if (weights_sum > 1 || weights_sum < 0)
    return;
  let zero_count = weights.reduce((ts, el)=>ts+=(el==0?1:0), 0);

  let remainder = zero_count > 0 ? (1 - weights_sum) / zero_count : 0;
  console.log(remainder);
  weights = weights.map(e=>e==0?remainder:e);

  prev_pfield = prev_pfield.map(e=>0);
  prev_cfield = prev_cfield.map(e=>0);
  let seed_index = 0;
  for (let i = 1; i <= n_cluster; i++) {
      desired_count[i] = pressure[i] = prev_pressure[i] = weights[i-1] * n_total  - 1;  // -1 because we have initialized a seed point
    //  seed_index = width * (i + 10) + 20 + 10*i;
      if (prev_pfield[seed_index] != 0) {
        while (!prev_pfield[(++seed_index) % n_total]) ;
      }
        //console.log("LOST " + i);
      prev_pfield[seed_index] = i;  // seed
      prev_cfield[seed_index] = 10;  // seed
      //seed_index +=  weights[i-1] * n_total; // Math.floor(Math.random() * Math.floor(width * height));
      seed_index =  Math.floor(Math.random() * Math.floor(width * height));
  }
  prev_pressure[0] = pressure[0] = -n_total + n_cluster;
  main();

}

var boost = n_cluster;

function rule(current, neighbors, ccurrent, cneighbors) {
  let own_sum = 0; let own_count = 0; let own_max = 0;
  let foreign_sum = 0; let foreign_count = 0; let foreign_max = 0;
  for (let i = 0; i < neighbors.length; i++)
    if (neighbors[i] == current) {
        own_sum += cneighbors[i];
        own_max = Math.max(own_max, cneighbors[i]);
        own_count++;
    } else {
        foreign_sum += cneighbors[i];
        foreign_max = Math.max(foreign_max, cneighbors[i]);
        foreign_count++;
    }

  let local_boost = (desired_count[current] - prev_pressure[current]) < 10 ? boost : 0;
  let new_ccurrent = current == 0 ? 0 :
//  (Math.max(prev_pressure[current], 0) / n_total + bias + local_boost + 0.5 * ccurrent + 0.5 * (own_sum - foreign_sum) / neighbors.length);// - 0.2 * foreign_sum / 4/*neighbors.length*/);// - (foreign_count ? foreign_sum / foreign_count : 0);
    prev_pressure[current] / n_total + bias + local_boost + 0.5 * ccurrent + 0.5 * (own_sum - foreign_sum) / neighbors.length;// - 0.2 * foreign_sum / 4/*neighbors.length*/);// - (foreign_count ? foreign_sum / foreign_count : 0);
  //if (new_ccurrent > n_total / n_cluster)
  //  new_ccurrent = n_total / n_cluster;
  if (new_ccurrent < 0)
    new_ccurrent =0 ;
  if (foreign_count == 0)
      return [current, new_ccurrent];

  let max_cn = AMaxI(cneighbors);
  let new_current = own_max >= foreign_max ? current : neighbors[max_cn[1]];
  if (new_current == 0) // Never let vacuum grow
    new_current = current;
  pressure[new_current]--;
  pressure[current]++;
  return [new_current, new_ccurrent];
}

function applyGraphRule(fr) {
    // Loop over all of the pixels
    for (var x=0; x<width; x++) {
        for (var y=0; y<height; y++) {
            // Get the pixel index
            let pixelindex = (y * width + x);
            // Constructing a 4-neighborhood graph
            let ne = [];
            let cne = [];
            if (x > 0) {
              ne.push(prev_pfield[pixelindex - 1]);
              cne.push(prev_cfield[pixelindex - 1]);
            }
            if (y > 0) {
              ne.push(prev_pfield[pixelindex - width]);
              cne.push(prev_cfield[pixelindex - width]);
            }
            if (x < width - 1) {
              ne.push(prev_pfield[pixelindex + 1]);
              cne.push(prev_cfield[pixelindex + 1]);
            }
            if (y < height - 1) {
              ne.push(prev_pfield[pixelindex + width]);
              cne.push(prev_cfield[pixelindex + width]);
            }
            let result = rule(prev_pfield[pixelindex], ne, prev_cfield[pixelindex], cne);
            pfield[pixelindex] = result[0];
            cfield[pixelindex] = result[1];
        }
    }

    for (let i = 0; i < cfield.length; i++) {
        prev_cfield[i] = (1 - urf_c) * prev_cfield[i] + urf_c * cfield[i];
        //if (prev_cfield[i] < 0) prev_cfield[i] = 0;
    }
    for (let i = 1; i < pressure.length; i++) {
        prev_pressure[i] = (1 - urf_p) * prev_pressure[i] + urf_p * pressure[i];
    }
    //prev_cfield.set(cfield);
    prev_pfield.set(pfield);
    if (pressure[0] == 0)
      bias = n_cluster / n_total;
    let ostr = "TS: " + fr + " Pressure: " + AMax(pressure)/*.join(", ") */+ " Max C: " + AMax(cfield);// + " Min C: " + AMin(cfield);
    document.getElementById("info").innerHTML = ostr;
}

// Create the image
function createImage(offset) {
    var cmax = AMax(cfield);
    var cmin = AMin(cfield);
    // Loop over all of the pixels
    for (var x=0; x<width; x++) {
        for (var y=0; y<height; y++) {
            // Get the pixel index
            var pixelindex = (y * width + x) * 4;

            // pfield
            var byteval = pfield[pixelindex / 4];// * 255 / (n_cluster + 1);

            // Set the pixel data
            imagedata1.data[pixelindex] = cb_p[byteval][0];     // Red
            imagedata1.data[pixelindex+1] = cb_p[byteval][1]; // Green
            imagedata1.data[pixelindex+2] = cb_p[byteval][2];  // Blue
            imagedata1.data[pixelindex+3] = 255;   // Alpha

            // cfield
            byteval = Math.floor(255 * (cfield[pixelindex / 4] - cmin) / (cmax - cmin));
            let cb = cb_c[byteval];
            if (!cb)
                console.log(byteval);
            // Set the pixel data
            imagedata2.data[pixelindex] = cb[0];     // Red
            imagedata2.data[pixelindex+1] = cb[1]; // Green
            imagedata2.data[pixelindex+2] = cb[2];  // Blue
            imagedata2.data[pixelindex+3] = 255;   // Alpha
        }
    }
}

ts = 0;
halt = false;
// Main loop
function main(tframe) {
    // Request animation frames
    //if (ts < 2000)
    if (!halt) window.requestAnimationFrame(main);
    for (let i = 0; i < 5; i++) {
        ts++;
        applyGraphRule(ts);
    }
    // Create the image
    createImage(ts);
    // Draw the image data to the canvas
    context1.putImageData(imagedata1, 0, 0);
    context2.putImageData(imagedata2, 0, 0);
}

// Call the main loop
restart();
