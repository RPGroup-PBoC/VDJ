
    // Determine the selection and assign to group filters.



    // Determine which point mutation is clicked
    var clicked_base = 0
    clicked_base = seq_source.data['mutant'][seq_source.selected['1d'].indices[0]];
    console.log(clicked_base)
    var mut = endog_sel.value;
    endog_filter.group = mut;
    target_filter.group = mut;

    // Update the sequence view
    seq_view.filters[0] = endog_filter;
    seq_source.data.view = seq_view;
    seq_source.change.emit();

    // Determine what point mutants should be plotted
    var points = []
    var colors = []
    for (var i = 0; i < seq_source.data['endogenous'].length ; i++) {
        var endo = seq_source.data['endogenous'][i]
        var mutation = seq_source.data['mutant'][i]
        if (endo === mut && mutation !== "No Mutation") {
            points.push(seq_source.data['mutant'][i]);
            colors.push(seq_source.data['color'][i]);
        }
    
    }

    // Define arrays for updating the endogenous plots
    var endog_data = [rep_endog, pooled_endog, dwell_all_endog, dwell_cut_endog, dwell_unloop_endog]
    var endog_views = [rep_endog_view, pooled_endog_view, dwell_all_endog_view, dwell_cut_endog_view, dwell_unloop_endog_view]

    // Push updates of the endogenous case
    for (var i = 0; i < endog_views.length; i++) { 
       endog_views[i].filters[0] = target_filter;
       endog_data[i].data.view = endog_views[i]
       endog_data[i].change.emit();}


    // Process the point data
    var point_data = [rep_point, pooled_point, dwell_unloop_point, dwell_cut_point, 
                      dwell_all_point];
    var point_views = [rep_point_view, pooled_point_view, dwell_unloop_point_view, 
                        dwell_cut_point_view, dwell_all_point_view];
    var point_filters = [rep_filter, pooled_filter, dwell_unloop_filter, 
                         dwell_cut_filter, dwell_all_filter];
    for (var i = 0; i < point_data.length; i++) { 
        var indices = [ ];
        var alphas = []
        for (var j = 0; j < point_data[i].data['mutant'].length; j++) {
            var point_mutant = point_data[i].data['mutant'][j];
            if (points.includes(point_mutant) === true) { 
                indices.push(j) 
                point_data[i].data['color'][j] = colors[points.indexOf(point_mutant)];
            }
            if (point_mutant === clicked_base | clicked_base === 0 | clicked_base === 'No Mutation' | clicked_base === undefined) {
                var set_alpha = 1
            }
            else { var set_alpha = 0.2}
            alphas.push(set_alpha)
        }
        // Add the indices to the corresponding filter. 
        point_filters[i].indices = indices;
        point_views[i].filters[0] = point_filters[i];
        point_data[i].data.view = point_views[i];
        point_data[i].data['alpha'] = alphas
        point_data[i].change.emit();
    }

    // Process dwell times
    var dwell_data = [dwell_all_point, dwell_cut_point, dwell_unloop_point];
    var dwell_views = [dwell_all_blank, dwell_cut_blank, dwell_unloop_blank]
    for (var i = 0; i < dwell_data.length; i++) { 
        var xs = [];
        var ys = [];
        var alphas = [];
        for (var j = 0; j < points.length; j++) {
            var target_point = points[j];
            var point_xs = [];
            var point_ys = [];
            if (target_point === clicked_base | clicked_base === 0 | clicked_base === 'No Mutation') {
                var set_alpha = 1;
            }
            else { var set_alpha=0.25;}
            for (var k = 0; k < dwell_data[i].data['mutant'].length; k++) { 
                var point_mutant = dwell_data[i].data['mutant'][k];

                    if (target_point === point_mutant) { 
                        point_xs.push(dwell_data[i].data['dwell_time'][k]);
                        point_ys.push(dwell_data[i].data['ecdf'][k]) ;
                    }

                   
            }
            xs.push(point_xs);
            ys.push(point_ys);
            alphas.push(set_alpha)
        }
        console.log(alphas)

 
        // Update the displayed data source with the correct multiline
        // parameters
        dwell_views[i].data['xs'] = xs;
        dwell_views[i].data['ys'] = ys;
        dwell_views[i].data['mutant'] = points;
        dwell_views[i].data['c'] = colors;
        dwell_views[i].data['alpha'] = alphas;
        dwell_views[i].change.emit();
    }
