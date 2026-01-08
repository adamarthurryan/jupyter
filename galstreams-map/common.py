import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Galactocentric
import astropy.units as u
from matplotlib.patches import Ellipse

# simplified galactic model
bulge_radius = 2 * u.kpc
disk_radius = 20 * u.kpc
disk_height = 0.5 * u.kpc 
stellar_halo_radius = 100 * u.kpc
dark_halo_radius = 292 * u.kpc


gc_frame = Galactocentric()
sun_coord=SkyCoord(ra=0*u.deg, dec=0*u.deg, distance=0*u.kpc).transform_to(gc_frame)

def plot_galaxy_map(ax_front, ax_top):

    # Disk and bulge patches
    disk_top = plt.Circle((0, 0), disk_radius.to(u.kpc).value, color='lightblue', alpha=0.35, label='Galactic Disk')
    #disk_front = plt.Rectangle((-disk_radius.to(u.kpc).value, -disk_height.to(u.kpc).value), 2*disk_radius.to(u.kpc).value, 2*disk_height.to(u.kpc).value, color='white', alpha=0.35, label='Galactic Disk')
    disk_front = Ellipse((0,0), 2*disk_radius.to(u.kpc).value, 2*disk_height.to(u.kpc).value, color='white', alpha=0.35, label='Galactic Disk')
    bulge_top = plt.Circle((0, 0), bulge_radius.to(u.kpc).value, color='white', alpha=0.65, label='Galactic Bulge')
    bulge_front = plt.Circle((0, 0), bulge_radius.to(u.kpc).value, color='white', alpha=0.65, label='Galactic Bulge')

    ax_front.add_patch(disk_front)
    ax_front.add_patch(bulge_front)
    ax_top.add_patch(disk_top)
    ax_top.add_patch(bulge_top)

    # sun
    ax_front.scatter(sun_coord.x.value, sun_coord.z.value, c='yellow', marker='*', s=50, alpha=1, zorder=999, label="Sun")
    ax_top.scatter(sun_coord.x.value, sun_coord.y.value, c='yellow', marker='*', s=50, alpha=1, zorder=999, label="Sun")

    # halos
    stellar_halo_top = plt.Circle((0, 0), stellar_halo_radius.to(u.kpc).value, color='darkblue', alpha=0.3, label='Stellar Halo')
    
    stellar_halo_front = plt.Circle((0, 0), stellar_halo_radius.to(u.kpc).value, color='darkblue', alpha=0.3, label='Stellar Halo')

    dark_halo_top = plt.Circle((0, 0), dark_halo_radius.to(u.kpc).value, alpha=1, facecolor='none', edgecolor='white', linestyle='--', label='Dark Halo')
    dark_halo_front = plt.Circle((0, 0), dark_halo_radius.to(u.kpc).value, alpha=1, facecolor='none', edgecolor='white', linestyle='--', label='Dark Halo')
    
    ax_front.add_patch(stellar_halo_front)
    ax_top.add_patch(stellar_halo_top)
    ax_front.add_patch(dark_halo_front)
    ax_top.add_patch(dark_halo_top)

# Set labels, title, and legend with light colors
def theme_galaxy_map(ax_front, ax_top):
    for ax in [ax_front, ax_top]:
        ax.set_facecolor('black')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.tick_params(axis='y', colors='white')
        #ax.legend(facecolor='black', edgecolor='white', labelcolor='white')
        ax.set_aspect('equal', 'box')
        ax.grid(True, color='gray', alpha=0.2)

    ax_front.set_xlabel('X (kpc)')
    ax_front.set_ylabel('Z (kpc)')
    ax_front.set_title('Edge-On (X-Z Plane)')
    ax_top.set_xlabel('X (kpc)')
    ax_top.set_ylabel('Y (kpc)')
    ax_top.set_title('Top-Down (X-Y Plane)')


def do_plot(streams_galactocentric, clusters_xyz, streams=False, clusters=False, mag_clouds=False, x_scale=25, yz_scale=15):
    clusters_x, clusters_y, clusters_z = clusters_xyz
    
    fig, [ax_front, ax_top] = plt.subplots(2, 1, figsize=(10,8), layout='constrained', sharex=True, facecolor='black')

    plot_galaxy_map(ax_front, ax_top)
    theme_galaxy_map(ax_front, ax_top)

    # Globular clusters
    if (clusters):
        ax_front.scatter(clusters_x, clusters_z, c='cyan', marker='o', s=10, label='Globular Clusters', zorder=10)
        ax_top.scatter(clusters_x, clusters_y, c='cyan', marker='o', s=10, label='Globular Clusters', zorder=10)

    if (streams):
        stream_colors = ['white']
        for idx, (name, stream_coords) in enumerate(streams_galactocentric.items()): 
            x_stream = stream_coords.x.value
            y_stream = stream_coords.y.value
            z_stream = stream_coords.z.value
            ax_front.plot(x_stream, z_stream,label="Stellar Stream", alpha=0.4, linestyle='-', color=stream_colors[idx % len(stream_colors)], linewidth=1.5)
            ax_top.plot(x_stream, y_stream, label="Stellar Stream", alpha=0.4, linestyle='-', color=stream_colors[idx % len(stream_colors)], linewidth=1.5)

    if (mag_clouds):
        # Magellanic Clouds
        lmc = SkyCoord(ra=80.89416667*u.deg, dec=-69.75611111*u.deg, distance=49.97*u.kpc).transform_to(gc_frame)
        smc = SkyCoord(ra=13.18666667*u.deg, dec=-72.82861111*u.deg, distance=62.44*u.kpc).transform_to(gc_frame)
        ax_front.scatter(lmc.x.value, lmc.z.value, c='orange', marker='o', s=100, label='LMC', zorder=10)
        ax_top.scatter(lmc.x.value, lmc.y.value, c='orange', marker='o', s=100, label='LMC', zorder=10)
        ax_front.scatter(smc.x.value, smc.z.value, c='red', marker='o', s=100, label='SMC', zorder=10)
        ax_top.scatter(smc.x.value, smc.y.value, c='red', marker='o', s=100, label='SMC', zorder=10)

    for ax in [ax_front, ax_top]:
        ax.set_xlim(-x_scale, x_scale)
        ax.set_ylim(-yz_scale, yz_scale)

    # create a shared legend for the plot
    handles1, labels1 = ax_top.get_legend_handles_labels()
    by_label = dict(zip(labels1, handles1))
    fig.legend(by_label.values(), by_label.keys(), loc='center right') # Adjust bbox_to_anchor for placement outside subplots

    plt.tight_layout()
    plt.show()
