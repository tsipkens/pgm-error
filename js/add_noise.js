// create an array with linear spacing
var linspace = function(a, b, n) {
  if (typeof n === "undefined") n = Math.max(Math.round(b - a) + 1, 1);
  if (n < 2) {
    return n === 1 ? [a] : [];
  }
  var i, ret = Array(n);
  n--;
  for (i = n; i >= 0; i--) {
    ret[i] = Math.round((i * b + (n - i) * a) / n * 100) / 100;
  }
  return ret;
};

// standard random normal numbers
var randn = function() {
  var u = 0,
    v = 0;
  while (u === 0) u = Math.random(); //Converting [0,1) to (0,1)
  while (v === 0) v = Math.random();
  return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

var add_noise = function(s_bar, tau, the, gam, N_shots, n) {

  N_s = 1 // s_bar.length // length of each signal

  s = s_bar

  for (ii = 0; ii < s_bar.length; ii++) {
    // Setup for random variables
    n_P = randn(); // standard normal random vector, realizes Poisson noise
    n_G = randn(); // standard normal random vector, realizes Gaussian noise

    // Generate observed signals by adding error terms
    s[ii] = s_bar[ii] * 1 + // expected average signal (s_bar)
      tau * s_bar[ii] * n + // shot-to-shot error (delta)
      Math.sqrt(the) * Math.sqrt(s_bar[ii] * 1 + tau * s_bar[ii] * n) * n_P + // Poisson noise (p)
      gam * n_G; // Gaussian noise (g)
  }

  return s
}



var N_shots = 1,
  gam_vec = linspace(0, 8, 21),
  the_vec = [0, 1 / 100, 1 / 20, 1 / 10, 1 / 7, 1 / 5, 1 / 4, 1 / 3, 1 / 2.5, 1 / 2,
    1 / 1.5, 1 / 1.2, 1, 1 / 0.7, 1 / 0.5, 1 / 0.4, 1 / 0.3, 1 / 0.2, 1 / 0.1
  ]
tau_vec = linspace(0, 0.4, 21);

// get parameters
var the = the_vec[document.getElementById('theSlider').value - 1],
  gam = gam_vec[document.getElementById('gamSlider').value - 1],
  tau = tau_vec[document.getElementById('tauSlider').value - 1];

// display current parameters
document.getElementById('theval').value = the;
document.getElementById('gamval').value = gam;
document.getElementById('tauval').value = tau;

// print to console
console.log([tau, gam, the])

var realize_noise = function(s, n) {
  return add_noise(s, tau, the, gam, N_shots, n)
}
var J_peak = function(the) {
  return Math.round(100 / the * 10) / 10;
}
document.getElementById('theval2').value = J_peak(the);

var n1 = -2,
  n2 = 0,
  n3 = 2;


// set the dimensions and margins of the graph
var margin = {
    top: 0,
    right: 60,
    bottom: 50,
    left: 60
  },
  width = 700 - margin.left - margin.right,
  height = 500 - margin.top - margin.bottom;

// append the svg object to the body of the page
var svg = d3.select("#my_dataviz")
  .append("svg")
  .attr("width", width + margin.left + margin.right)
  .attr("height", height + margin.top + margin.bottom)
  .append("g")
  .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

d3.csv("https://raw.githubusercontent.com/tsipkens/wat-lii-error/master/data/gaus.csv", function(data) {

  // Add X axis
  var x = d3.scaleLinear()
    .domain([0, 2500])
    .range([0, width]);
  svg.append("g")
    .attr("transform", "translate(0," + height + ")")
    .call(d3.axisBottom(x));
  svg.append("g")
    .call(d3.axisTop(x));

  // Add Y axis
  var y = d3.scaleLinear()
    .domain([-0.30, 1.95])
    .range([height, 0]);
  svg.append("g")
    .call(d3.axisLeft(y));
  var yvr = d3.scaleLinear()
    .domain([100 / the * (-0.30), 100 / the * 1.95])
    .range([height, 0]);
  var yAxis2 = svg.append("g")
    .attr("transform", "translate(" + width + ",0)")
    .call(d3.axisRight(yvr))

  //-- Add axis labels --//
  // Add X axis label:
  svg.append("text")
    .attr("text-anchor", "middle")
    .attr('x', width / 2)
    .attr('y', height + 35)
    .text("Time, t");

  // Y axis label:
  svg.append("text")
    .attr("text-anchor", "middle")
    .attr('transform', 'translate(-35,' + height / 2 + ')rotate(-90)')
    .text("Scaled signal, s/max{s}")
  svg.append("text")
    .attr("text-anchor", "middle")
    .attr('transform', 'translate(' + (52 + width) + ',' + height / 2 + ')rotate(-90)')
    .text("Signal, s [counts]")

  // Fill in the main plot ---------------------------------------------------//
  svg.append("path")
    .datum(data)
    .attr("fill", "none")
    .attr("stroke", "#E40066")
    .attr("stroke-width", 2)
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.t)
      })
      .y(function(d) {
        return y([d.J] / 100)
      })
    )

  svg.append("path")
    .datum(data)
    .attr("fill", "none")
    .attr("stroke", "#EAC435")
    .attr("stroke-width", 1.5)
    .attr('stroke-dasharray', "4 3")
    .attr("id", 'l1')
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.t)
      })
      .y(function(d) {
        return y((1 + tau * n1) * [d.J] / 100)
      })
    )

  svg.append("path")
    .datum(data)
    .attr("fill", "none")
    .attr("stroke", "#0E1C36")
    .attr("stroke-width", 1.5)
    .attr('stroke-dasharray', "4 3")
    .attr("id", 'l2')
    .attr("d", d3.line()
      .x(function(d) {
        return x(d.t)
      })
      .y(function(d) {
        return y((1 - tau * n1) * [d.J] / 100)
      })
    )

  svg.append('g')
    .selectAll("dot")
    .data(data)
    .enter()
    .append("circle")
    .attr("cx", function(d) {
      return x(d.t);
    })
    .attr("cy", function(d) {
      return y(realize_noise([d.J], n1) / 100);
    })
    .attr("r", 2)
    .attr("stroke", "black")
    .attr("stroke-width", 0)
    .attr("id", 'A')
    .style("fill", '#EAC435')

  svg.append('g')
    .selectAll("dot")
    .data(data)
    .enter()
    .append("circle")
    .attr("cx", function(d) {
      return x(d.t);
    })
    .attr("cy", function(d) {
      return y(realize_noise([d.J], n2) / 100);
    })
    .attr("r", 2)
    .attr("stroke", "black")
    .attr("stroke-width", 0)
    .attr("id", 'B')
    .style("fill", '#E40066')

  svg.append('g')
    .selectAll("dot")
    .data(data)
    .enter()
    .append("circle")
    .attr("cx", function(d) {
      return x(d.t);
    })
    .attr("cy", function(d) {
      return y(realize_noise([d.J], n3) / 100);
    })
    .attr("r", 2)
    .attr("stroke", "black")
    .attr("stroke-width", 0)
    .attr("id", 'C')
    .style("fill", '#0E1C36')

  // slider controls
  d3.select("#gamSlider").on("change", function() {
    gam = gam_vec[this.value - 1]
    updatePlot()
  })
  d3.select("#theSlider").on("change", function() {
    the = the_vec[this.value - 1]
    updatePlot()
  })
  d3.select("#tauSlider").on("change", function() {
    tau = tau_vec[this.value - 1]
    updatePlot()
  })

  //------------------------------------------------------------------------//
  // a generic plot updater using a given data set
  // e.g., used whenever the slider for da is changed
  function updatePlot() {
    var realize_noise = function(s, n) {
      return add_noise(s, tau, the, gam, N_shots, n)
    }

    // give these new data to update plot
    svg.select("#l1")
      .datum(data)
      .transition()
      .duration(100)
      .attr("d", d3.line()
        .x(function(d) {
          return x(d.t)
        })
        .y(function(d) {
          return y((1 + tau * n1) * [d.J] / 100)
        })
      );

    svg.select("#l2")
      .datum(data)
      .transition()
      .duration(100)
      .attr("d", d3.line()
        .x(function(d) {
          return x(d.t)
        })
        .y(function(d) {
          return y((1 - tau * n1) * [d.J] / 100)
        })
      );

    svg.selectAll("#A")
      .data(data)
      .transition()
      .duration(100)
      .attr("cy", function(d) {
        return y(realize_noise([d.J], n1) / 100);
      });

    svg.selectAll("#B")
      .data(data)
      .transition()
      .duration(100)
      .attr("cy", function(d) {
        return y(realize_noise([d.J], n2) / 100);
      });

    svg.selectAll("#C")
      .data(data)
      .transition()
      .duration(100)
      .attr("cy", function(d) {
        return y(realize_noise([d.J], n3) / 100);
      });

    var yvr = d3.scaleLinear()
      .domain([100 / the * (-0.30), 100 / the * 1.95])
      .range([height, 0]);
    yAxis2.attr("transform", "translate(" + width + ",0)")
      .call(d3.axisRight(yvr))
  }
  //------------------------------------------------------------------------//

  // setpoint to show Poisson noise variations due to shot-to-shot
  d3.select("#set2").on("click", function() {
    setvals(1, 19, 8);
  })

  d3.select("#set1").on("click", function() {
    setvals(1, 12, 21);
  })

  // function to set sliders to specific values
  function setvals (gamv, thev, tauv) {
    gam = gam_vec[gamv - 1]; the = the_vec[thev - 1]; tau = tau_vec[tauv - 1];

    document.getElementById('gamval').value = Math.round(gam * 100) / 100;
    document.getElementById('theval').value = Math.round(the * 100) / 100;
    document.getElementById('tauval').value = Math.round(tau * 100) / 100;
    document.getElementById("theval2").value = J_peak(tau);

    d3.select("#gamSlider")
      .transition().attr("value", gamv)
    d3.select("#theSlider")
      .transition().attr("value", thev)
    d3.select("#tauSlider")
      .transition().attr("value", tauv)

    updatePlot()
  }

})


function displayval (val, vec, id) {
  document.getElementById(id).value = Math.round(vec[val - 1] * 100) / 100;
  if (id == "theval") {
    document.getElementById("theval2").value = J_peak(vec[val - 1]);
  }
}
