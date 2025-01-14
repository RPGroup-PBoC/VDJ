import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import vdj.viz
import vdj.io
vdj.viz.plotting_style()

# Load all cutting probability estimates taking gaussian approximation.
cut_data = pd.read_csv('../../data/pooled_cutting_probability.csv', comment='#')
cut_data = cut_data[(cut_data['hmgb1'] == 80) & (cut_data['salt']=='Mg') & (cut_data['mutant']=='WT12rss')]

# Load the precomputed posterior distributioons
cut_posts = pd.read_csv('../../data/pooled_cutting_probability_posteriors.csv', 
                        comment='#')
cut_posts = cut_posts[(cut_posts['hmgb1']==80) & (cut_posts['salt']=='Mg') & (cut_posts['mutant']=='WT12rss')]

#%%
fig, ax = plt.subplots(1, 1, figsize=(4.1, 2.3))
ax.fill_between(cut_posts['probability'], 0, cut_posts['posterior'],
                color='slategrey', alpha=0.4)
ax.plot(cut_posts['probability'], cut_posts['posterior'], color='white')
ax.plot(cut_data['mean'], cut_posts['posterior'].max()/3, zorder=10,
            color='dodgerblue', marker='o', ms=15, markerfacecolor='white')
ax.hlines(cut_posts['posterior'].max()/3, cut_data['mean'] - cut_data['std'],
        cut_data['mean'] + cut_data['std'], color='dodgerblue', zorder=10,
        lw=4)

_ = ax.set_xlim([0.3, 0.7])
_ = ax.set_ylim([0.0, cut_posts['posterior'].max()+0.005])
_ = ax.set_ylabel('$P(p_\mathrm{cut} | n_\mathrm{loops}, n_\mathrm{cuts})$',
                fontsize=15)
_ = ax.set_xlabel(r'$p_\mathrm{cut}$', fontsize=15)

bead_stats = '\n'.join((
    'V4-57-1 12RSS',
    '(reference)',
    r'$n_\mathrm{loops}=%i$' % (cut_data['n_loops'], ),
    r'$n_\mathrm{cuts}=%i$' % (cut_data['n_cuts'], )))

bbox_props = dict(boxstyle='square', edgecolor='k', facecolor='white', alpha=0.5)

#_ = ax.text(0.66, 0.051, bead_stats, fontsize=12, horizontalalignment='center',
#            verticalalignment='top', bbox=bbox_props)
_ = ax.text(cut_data['mean'], cut_posts['posterior'].max()/3 + 0.008,
            r'$\mu \pm \sigma$', fontsize=12, verticalalignment='center',
            horizontalalignment='center')
_ = ax.text(cut_data['mean'], cut_posts['posterior'].max()/3-0.0005, 'N',
            fontsize=13, color='dodgerblue', verticalalignment='center', 
            horizontalalignment='center', zorder=15)

fig.savefig('./SubFigX_point_posterior_definition.pdf', bbox_inches='tight',
            facecolor='white')

# %%
