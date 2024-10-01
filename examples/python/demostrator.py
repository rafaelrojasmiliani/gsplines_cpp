
import matplotlib.pyplot as plt
import numpy as np
import copy
from matplotlib.widgets import RadioButtons, Button, Slider
import inspect

try:
    import gsplines
    import gsplines.plot
    import gsplines.optimization
except ImportError:
    import pathlib
    import sys
    cwd = pathlib.Path(__file__).parent.absolute()
    setup_script = pathlib.Path(
        cwd, '../../', 'setup_python_env.py')

    if not setup_script.exists():
        print('Import Error: Seems that you haven\'t compile the project')
        sys.exit(1)

    exec(setup_script.read_text())
    try:
        import gsplines
        import gsplines.plot
        import gsplines.optimization

    except ImportError:
        print('Import Error: Seems that you failed to compile properly')
        sys.exit(1)


def main():
    global waypoints
    global prev_waypoints
    global optimization_function

    waypoints = []
    prev_waypoints = None
    trj_array = []
    slider = None
    optimization_function_array =\
        [gsplines.optimization.broken_lines_path,
         gsplines.optimization.minimum_acceleration_path,
         gsplines.optimization.minimum_jerk_path,
         gsplines.optimization.minimum_snap_path,
         gsplines.optimization.minimum_crackle_path,
         lambda x: gsplines.optimization.rojas_path(x, slider.val)]
    optimization_function = optimization_function_array[0]

    fig = plt.figure(figsize=(10, 6))

    ax_big = plt.subplot2grid((4, 5), (0, 0), rowspan=4, colspan=3)
    ax_big.set_xlim(-1, 1)
    ax_big.set_ylim(-1, 1)

    axes_small = np.array([[plt.subplot2grid((4, 5), (i, j+3))
                            for j in range(2)] for i in range(4)])

    def on_click(event):
        toolbar = plt.get_current_fig_manager().toolbar
        # Only respond to clicks in the big axis
        # print(toolbar.mode)
        if event.inaxes == ax_big \
                and toolbar.mode not in ['zoom rect', 'pan/zoom']:
            waypoints.append([event.xdata, event.ydata])
            wp = np.array(waypoints)
            ax_big.plot(wp[:, 0],
                        wp[:, 1], 'ro')  # 'o' marker for circles
            ax_big.set_xlim(-1, 1)
            ax_big.set_ylim(-1, 1)
            fig.canvas.draw_idle()
            fig.canvas.blit(ax_big.bbox)
            fig.canvas.flush_events()

    fig.canvas.mpl_connect('button_press_event', on_click)

# Radio buttons on the big axis
    # Position for radio buttons

# Function to handle radio button selection
    def on_radio_change(label):
        radio_id = radio_button_labels.index(label)
        global optimization_function
        optimization_function = optimization_function_array[radio_id]

    radio_button_labels = [
        'min speed', 'min acceleration',
        'min jerk', 'min snap',
        'min crackle (unstable)',
        'speed jerk balance']

    # Position for the slider
    ax_slider = plt.axes([0.005, 0.1, 0.01, 0.8], facecolor="lightgray")
    slider = Slider(ax_slider, 'Range', valmin=0,
                    valmax=10, valinit=1.5, orientation='vertical')

    ax_radio = plt.axes([0.02, 0.65, 0.1, 0.2],
                        facecolor="lightgoldenrodyellow")
    ax_button = plt.axes([0.02, 0.55, 0.1, 0.05])  # Position for the button
    ax_clear_button = plt.axes([0.02, 0.45, 0.1, 0.05])
    ax_reuse_button = plt.axes([0.02, 0.35, 0.1, 0.05])

    radio = RadioButtons(ax_radio, radio_button_labels)
    radio.on_clicked(on_radio_change)
    button = Button(ax_button, 'Calculate')
    clear_button = Button(ax_clear_button, 'Clear')
    reuse_button = Button(ax_reuse_button, 'Use previous Waypoins')

# Function to handle button click

    def on_button_click(event):
        global prev_waypoints
        global waypoints
        wp = np.array(waypoints)
        if wp.shape[0] < 2:
            print("cannto compute trajectory less than two points")
            return
        if not optimization_function:
            print('optimization function is emtpy')
            return

        gspline = optimization_function(wp)

        if not gspline:
            prev_waypoints = copy.deepcopy(waypoints)
            waypoints.clear()
            print('Optimization process Failed!')
        else:
            [ax.cla() for ax in [ax_big] + axes_small.flatten().tolist()]
            wp = np.array(waypoints)
            trj_array.append(gspline)
            for i, trj in enumerate(trj_array):
                time_spam = np.arange(*trj.get_domain(), 0.001)
                points_to_plot = trj(time_spam)
                ax_big.plot(points_to_plot[:, 0], points_to_plot[:, 1])
                ax_big.set_xlim(-1, 1)
                ax_big.set_ylim(-1, 1)

            for i, trj in enumerate(trj_array):
                gsplines.plot.plot_derivatives_in_axes(
                    trj, axes_small, color=ax_big.get_lines()[i].get_color())
            ax_big.plot(wp[:, 0],
                        wp[:, 1], 'ro')  # 'o' marker for circles
            prev_waypoints = copy.deepcopy(waypoints)
            waypoints.clear()
        fig.canvas.draw_idle()
        fig.canvas.blit(ax_big.bbox)
        fig.canvas.flush_events()

    def on_clear_click(event):
        ax_big.cla()  # Clear the big axis
        ax_big.set_title("Big Axis (Click Here)")  # Restore the title
        [ax.cla() for ax in [ax_big] + axes_small.flatten().tolist()]
        ax_big.set_xlim(-1, 1)
        ax_big.set_ylim(-1, 1)
        trj_array.clear()
        waypoints.clear()
        fig.canvas.draw_idle()
        fig.canvas.flush_events()

    def on_reuse_buttom(event):
        global waypoints
        global prev_waypoints
        print(prev_waypoints)
        waypoints = copy.deepcopy(prev_waypoints)
        on_button_click(event)

    button.on_clicked(on_button_click)
    clear_button.on_clicked(on_clear_click)
    reuse_button.on_clicked(on_reuse_buttom)
    plt.show()


if __name__ == "__main__":
    main()
# # Event handler for mouse clicks on the big axis only


# # Connect the click event to the figure

# # Add some labels for visualization purposes
# ax_big.set_title("Big Axis (Click Here)")
# for i in range(2):
#     for j in range(4):
#         axes_small[i, j].set_title(f"Small Axis {i*4 + j + 1}")

# # Radio buttons on the big axis
# # Position for radio buttons
# ax_radio = plt.axes([0.05, 0.65, 0.1, 0.2], facecolor="lightgoldenrodyellow")
# radio = RadioButtons(ax_radio, [
#                      'Option 1', 'Option 2', 'Option 3', 'Option 4', 'Option 5', 'Option 6'])

# # Function to handle radio button selection


# def on_radio_change(label):
#     print(f"Radio button selected: {label}")


# # Connect radio button event
# radio.on_clicked(on_radio_change)


# # Vertical RangeSlider on the big axis
# # Position for the slider
# ax_slider = plt.axes([0.15, 0.1, 0.03, 0.8], facecolor="lightgray")
# slider = RangeSlider(ax_slider, 'Range', valmin=0, valmax=10,
#                      valinit=(2, 8), orientation='vertical')

# # Function to handle slider change


# def on_slider_change(val):
#     print(f"Range selected: {val}")


# # Connect the slider change event
# slider.on_changed(on_slider_change)

# # Display the plot
# plt.tight_layout()
