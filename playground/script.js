// Helper
Float64Array.prototype.max = Array.prototype.max = function() {
  return this.reduce((a, b)=> Math.max(a, b));
}

Float64Array.prototype.min = Array.prototype.min = function() {
  return this.reduce((a, b)=> Math.min(a, b));
}

Float64Array.prototype.sum = Array.prototype.sum = function() {
  return this.reduce((a, b)=> a+b);
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
var bias = 0.1;//n_cluster / n_total;//1/(width*height);

var pfield = new Uint32Array(width * height);
var cfield = new Float64Array(width * height);
var pressure = new Array(n_cluster + 1);
var desired_count = new Array(n_cluster + 1);
var endangered = new Array(n_cluster + 1);

function sin_col(i, size, phase) {
  var sin = Math.sin(Math.PI / (size / 2) * i + phase);
  var int = Math.floor(sin * 127) + 128;
  return int;
}

var cb_c = new Array(256);
for (var i = 0; i < 256; i++) {
  var red   = Math.min(i/100,1) * sin_col(i, 255, 0 * Math.PI * 2/3); // 0   deg
  var blue  = Math.min(i/100,1) * sin_col(i, 255, 1 * Math.PI * 2/3); // 120 deg
  var green = Math.min(i/100,1) * sin_col(i, 255, 2 * Math.PI * 2/3); // 240 deg
  cb_c[i] = [red, green, blue];
}

// Init
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
 cb_p[1] = [155,155,155]; 
 cb_p[2] = [128,128,128]; 


  w0 = Number(document.getElementById("w0").value)
  w1 = Number(document.getElementById("w1").value)
  pressure = new Array(n_cluster + 1);
  weights = new Array(n_cluster);
  weights.fill(0);
  weights[0] = w0;
  weights[1] = w1;
  // auto balance
  let weights_sum = weights.sum()
  if (weights_sum > 1 || weights_sum < 0)
    return;
  let zero_count = weights.reduce((ts, el)=>ts+=(el==0?1:0), 0);

  let remainder = zero_count > 0 ? (1 - weights_sum) / zero_count : 0;
  console.log(remainder);
  weights = weights.map(e=>e==0?remainder:e);

  bias = 1/(weights.min()*n_total);

  pfield.fill(0);
  cfield.fill(0);

  let seed_index = 0;
  pressure[0] = 0; desired_count[0] = 0;
  for (let i = 1; i <= n_cluster; i++) {
      seed_index =  Math.floor(Math.random() * Math.floor(width * height));
      while (pfield[seed_index] != 0)
        seed_index =  Math.floor(Math.random() * Math.floor(width * height));

      pfield[seed_index] = i;  // seed
      cfield[seed_index] = 10;  // seed

      desired_count[i] = pressure[i] = Math.round(weights[i-1] * n_total) ;  // -1 because we have initialized a seed point
      pressure[i]--;
      //seed_index += Math.floor(desired_count[i]) ; // Math.floor(Math.random() * Math.floor(width * height));
  }
  pressure[1] += n_total - desired_count.sum();
  desired_count[1] += n_total - desired_count.sum();

  pressure[0] = -n_total + n_cluster;
  endangered.fill(true);endangered[0] = false;
  halt = true;
  setTimeout(main, 1000);
  setTimeout(()=>halt=false,900);

}

function applyGraphRule(fr) {
    // Loop over all pixels
    for (let pixelindex = 0; pixelindex < n_total; pixelindex++) {
        let ne = graph[pixelindex];
        let current = pfield[pixelindex];
        let ccurrent = cfield[pixelindex];

        let own_sum = 0; let own_count = 0; let own_max = 0; 
        let foreign_sum = 0; let foreign_count = 0; let foreign_max = 0; let foreign_max_ind = -1;
        for (let ne_idx of ne) {
          let cne = cfield[ne_idx];
          if (pfield[ne_idx] == current) {
              own_sum += cne;
              own_max = Math.max(own_max, cne);
              own_count++;
          } else {
              foreign_sum += cne;
              if (cne > foreign_max) {
                foreign_max = cne;
                foreign_max_index = ne_idx;
              }
              foreign_count++;
          }
        }

        let power = current == 0 ? 0 : (pressure[current] / n_total) / (weights[current-1]) + bias;
        if (power < 0)
          power = 0;
        let new_ccurrent = current == 0 ? 0 :
          power + 0.5 * ccurrent + 0.5 * (own_sum - foreign_sum) / ne.length;// - 0.2 * foreign_sum / 4/*neighbors.length*/);// - (foreign_count ? foreign_sum / foreign_count : 0);
        if (new_ccurrent < 0)
          new_ccurrent = 0;        
        cfield[pixelindex] = new_ccurrent;
        if (foreign_count == 0 || endangered[current]) { // Shortcut
          continue;
        }

        let new_current = own_max >= foreign_max ? current : pfield[foreign_max_index];
        if (new_current == 0) // Never let vacuum grow
          continue;
        
        pressure[new_current]--;
        pressure[current]++;
        pfield[pixelindex] = new_current;
    }

    for (let i = 1; i <= n_cluster; i++)
      endangered[i] = ((desired_count[i] - pressure[i]) < (0.5 * desired_count[i]));
}

// Create the image
function createImage(offset) {
    var cmax = cfield.max();
    var cmin = cfield.min();
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
halt = false ;

// Loop over all pixels to create graph
graph = [];
console.log(width + " " + height);
for (var y=0; y<height; y++) {
  for (var x=0; x<width; x++) {
        // Get the pixel index
        let pixelindex = (y * width + x);
        // Constructing a 4-neighborhood graph
        let ne = [];
        if (x > 0) {
          ne.push(pixelindex - 1);
        }
        if (y > 0) {
          ne.push(pixelindex - width);
        }
        if (x < width - 1) {
          ne.push(pixelindex + 1);
        }
        if (y < height - 1) {
          ne.push(pixelindex + width);
        }
        graph.push(ne);
    }
}
// Main loop
function main(tframe) {

    // Request animation frames
    //if (ts < 2000)
    if (!halt) 
      window.requestAnimationFrame(main);
//    for (let i = 0; i < (tframe ? tframe : 5); i++) {
    for (let i = 0; i < 50; i++) {
      ts++;
      applyGraphRule(ts);
    }
    let ostr = "Iter: " + ts + " Max imbalance (nodes): " + pressure.max().toFixed(2)/*.join(", ") */+ " Max Potential: " + cfield.max().toFixed(2) + " Avg Potential: " + cfield.sum() / n_total;// + " Min C: " + AMin(cfield);
    document.getElementById("info").innerHTML = ostr;
    //
    //for (let i = 1; i <= n_cluster; i++)
    //  if (desired_count[i] - prev_pressure[i] < 1 )
    //    restart();//alert(ts + ":  Partition " + i + " died.");
    // Create the image
    createImage(ts);
    // Draw the image data to the canvas
    context1.putImageData(imagedata1, 0, 0);
    context2.putImageData(imagedata2, 0, 0);
}

// Call the main loop
restart();
