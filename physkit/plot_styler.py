# Physkit - A Python toolkit for constants, unit conversions, and equations.
# Copyright (C) 2024 sapphimars
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.


import matplotlib.pyplot as plt


def plot_styler(
    xdata,
    ydata,
    title: str = None,
    xlabel: str = None,
    ylabel: str = None,
    xlims: list = None,
    ylims: list = None,
    color: str = None,
    label: str = None,
    marker: str = None,
    line: str = None,
    markersize: int = 4,
    linewidth: float = None,
    yerr: list = None,
    minorgrid: bool = False,
    loglog: bool = False,
    ax=None,
    scatter: bool = False,
) -> None:
    """
    Helper function which adds styling to each plot in the way we want.

    Args:
        xdata (list or np.ndarray),
        ydata (list or np.ndarray),
        title (str optional),
        xlabel (str optional),
        ylabel (str optional),
        xlims (list optional),
        ylims (list optional),
        color (str optional),
        label (str optional),
        marker (str optional),
        line (str optional)
        markersize (int default = 4),
        linewidth (float optional),
        yerr (np.ndarray optional),
        minorgrid (bool default = False),
        loglog (bool default = False),
        ax (matplotlib.axes.Axes optional),
        scatter (bool default = False),
    """

    if ax is None:
        ax = plt  # Default to using the global plt object if no axis

    if loglog:
        ax.set_xscale("log")
        ax.set_yscale("log")

    if yerr is not None:
        ax.errorbar(
            xdata,
            ydata,
            fmt=marker,
            yerr=yerr,
            color=color,
            markersize=markersize,
            label=label,
            linewidth=linewidth,
        )
    elif scatter:
        ax.scatter(xdata, ydata, marker=marker, color=color, label=label)
    else:
        ax.plot(
            xdata,
            ydata,
            marker=marker,
            color=color,
            markersize=markersize,
            label=label,
            linewidth=linewidth,
        )

    if title:
        ax.set_title(title, fontsize=15)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=12)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12)
    if xlims:
        ax.set_xlim(xlims)
    if ylims:
        ax.set_ylim(ylims)

    ax.grid(which="major", color="#CCCCCC")
    if minorgrid:
        ax.minorticks_on()
        ax.grid(which="minor", color="#EEEEEE", linestyle=":")

    ax.set_axisbelow(True)

    if label:
        ax.legend()
