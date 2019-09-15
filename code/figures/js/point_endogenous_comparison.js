
    // Determine which point mutation is clicked
    var clicked_base = 0
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
            if (hover_mut !== 'None' & hover_mut!== 'No Mutation' & hover_mut !== undefined) {
                if (mutation === hover_mut) { 
                    console.log('hover mut found')
                    colors.push(seq_source.data['color'][i])
                }
                else { 
                    console.log('hover mut not found')
                colors.push('#c2c2c2');
                }
            }
        else { 
            console.log('hover mut not selected')
            colors.push(seq_source.data['color'][i])
            }
        }
    
    }

    // Define arrays for updating the endogenous plots
    var endog_data = [rep_endog, pooled_endog, dwell_all_endog, 
                      dwell_cut_endog, dwell_unloop_endog, post_endog]
    var endog_views = [rep_endog_view, pooled_endog_view, dwell_all_endog_view, 
                        dwell_cut_endog_view, dwell_unloop_endog_view,
                        post_endog_view]

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
                point_data[i].data['display_color'][j] = colors[points.indexOf(point_mutant)];
            }
            if (point_mutant === clicked_base | clicked_base === 0 | clicked_base === 'No Mutation' | clicked_base === undefined) {
                var set_alpha = 1
            }
            else { var set_alpha = 1}
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
    var multi_data = [dwell_all_point, dwell_cut_point, dwell_unloop_point, post_point];
    var multi_views = [dwell_all_blank, dwell_cut_blank, dwell_unloop_blank, post_blank]
    for (var i = 0; i < multi_data.length; i++) { 
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
            else { var set_alpha=1;}
            for (var k = 0; k < multi_data[i].data['mutant'].length; k++) { 
                var point_mutant = multi_data[i].data['mutant'][k];

                    if (target_point === point_mutant) { 
                        point_xs.push(multi_data[i].data['x_val'][k]);
                        point_ys.push(multi_data[i].data['y_val'][k]) ;
                    }

                   
            }
            xs.push(point_xs);
            ys.push(point_ys);
            alphas.push(set_alpha)
        }

 
        // Update the displayed data source with the correct multiline
        // parameters
        multi_views[i].data['xs'] = xs;
        multi_views[i].data['ys'] = ys;
        multi_views[i].data['mutant'] = points;
        multi_views[i].data['c'] = colors;
        multi_views[i].data['alpha'] = alphas;
        multi_views[i].change.emit();

    }
