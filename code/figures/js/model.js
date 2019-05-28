// Variable definition
var prob_source = prob_source.data;
var dwell_time = dwell_source.data['time'];
var dwell_cdf = dwell_source.data['cdf'];
var dwell_pdf = dwell_source.data['pdf'];
var r_cut = rcut_slider.value;
var k_loop = kloop_slider.value;
var k_unloop = kunloop_slider.value;
var r = r_cut + k_unloop;


// Compute the quantities of interest
prob_source['p_loop'] = k_loop / (k_loop + k_unloop);
prob_source['p_cut'] = r_cut / (r_cut + k_unloop);
prob_source['f_loop'] = k_loop / (r + k_loop); 

// Compute the dwell time cumulative distribution
for (var i = 0; i < dwell_cdf.length; i++) {
    dwell_pdf[i] = r * Math.exp(-r * dwell_time[i]) 
    dwell_cdf[i] = 1 - Math.exp(-r * dwell_time[i]);
} 
dwell_pdf = sum(dwell)
prob_source.change.emit();
dwell_source.change.emit();