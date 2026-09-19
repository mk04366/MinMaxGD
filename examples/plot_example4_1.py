"""Plot references and outputs for Example 4.1.

The paper's reference/output comparison is Figure 4.3; Figure 4.4 is the
shared-resource usage schedule. Counter terms g^n d^t are represented here
as n firings completed by time t.

Run from the repository root:
    python examples/plot_example4_1.py

The script writes examples/example4_1_reference_vs_output.png.
"""

import csv
from pathlib import Path

import matplotlib.pyplot as plt



def collect_counter_data(csv_path):
    """Read the counters exported by the most recent C++ example run."""
    data = {subsystem: {"reference": [], "output": []}
            for subsystem in ("S1", "S2", "S3")}
    with csv_path.open(newline="") as file:
        for row in csv.DictReader(file):
            data[row["subsystem"]][row["kind"]].append(
                (int(row["gamma"]), None if row["delta"] == "inf" else int(row["delta"]))
            )
    return data


def counter_steps(events, horizon, plot_end):
    """Convert finite counter terms into their gamma-versus-delta staircase."""
    finite_events = [event for event in events if event[1] is not None]
    infinity_events = [event for event in events if event[1] is None]
    times = [0]
    counts = [0]

    for gamma, time in finite_events:
        times.append(time)
        counts.append(gamma)

    # The infinity term stores the final firing count beyond the finite horizon.
    if infinity_events:
        times.append(horizon)
        counts.append(infinity_events[0][0])

    # Extend the last known value to the derived plot end as visual padding.
    if times[-1] < horizon:
        times.append(horizon)
        counts.append(counts[-1])
    times.append(plot_end)
    counts.append(counts[-1])
    return times, counts


def plot_subsystem(axis, name, data, horizon, plot_end):
    """Plot one subsystem's reference and computed output schedules."""
    reference_times, reference_counts = counter_steps(data[name]["reference"], horizon, plot_end)
    output_times, output_counts = counter_steps(data[name]["output"], horizon, plot_end)

    axis.step(
        reference_times,
        reference_counts,
        where="post",
        color="#c0392b",
        linewidth=2.2,
        label="Reference",
    )
    axis.step(
        output_times,
        output_counts,
        where="post",
        color="#1769aa",
        linewidth=2.2,
        label="Output",
    )

    # Mark the finite counter transitions.
    axis.scatter(reference_times[1:-1], reference_counts[1:-1], color="#c0392b", s=28)
    # axis.scatter(output_times[1:-1], output_counts[1:-1], color="#1769aa", s=28)
    axis.set_title(f"Subsystem {name[1]}: reference vs output")
    axis.set_ylabel("Completed firings")
    axis.set_xlim(0, plot_end)
    axis.set_ylim(bottom=0)
    axis.grid(True, alpha=0.25)
    axis.legend(loc="upper left")


def main():
    """Create and save the three-subsystem comparison figure."""
    csv_path = Path(__file__).parents[1] / "misc" / "example4_1_plot_data.csv"
    if not csv_path.exists():
        raise FileNotFoundError(
            "Run build\\Example4_1.exe first; it exports the current counter data."
        )
    data = collect_counter_data(csv_path)
    finite_times = [
        time
        for subsystem in data.values()
        for series in subsystem.values()
        for _, time in series
        if time is not None
    ]
    horizon = max(finite_times, default=0)
    padding = max(5, min(10, max(1, horizon // 10)))
    plot_end = horizon + padding
    output_path = Path(__file__).parents[1] / "plots" / "example4_1_reference_vs_output.png"
    figure, axes = plt.subplots(3, 1, figsize=(11, 9), sharex=True)

    for axis, subsystem in zip(axes, data):
        plot_subsystem(axis, subsystem, data, horizon, plot_end)

    axes[-1].set_xlabel("Time")
    figure.suptitle("Example 4.1: reference counters and computed outputs", fontsize=15)
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(output_path, dpi=160)
    print(f"Saved plot to {output_path}")


if __name__ == "__main__":
    main()
