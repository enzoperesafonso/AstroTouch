"""Help shown in the GUI: one entry per settings section, plus the getting-started guide and tool hints."""

GETTING_STARTED = """
<h2>Getting started</h2>
<p>AstroTouch turns an astronomical image into a plaque you can feel: <b>height follows brightness</b>,
with a braille or raised-print label.</p>
<ol>
<li><b>Open</b> a FITS file or a picture (PNG, JPEG), or drag one onto the window.</li>
<li><b>Pick templates</b> at the top of the settings: a <i>plaque format</i> (size, border, label strip) and an
<i>image recipe</i> (how the image is cleaned and stretched). They are good starting points; you can change
anything afterwards.</li>
<li><b>Look at the 2D view.</b> It shows a shaded preview of the model and updates a moment after any change.
The chart under <i>Image and stretch</i> shows what the stretch does to your image.</li>
<li><b>Fix things by hand</b> with the tools above the 2D view: paint out anything unwanted, protect what must
stay, add textures, markers, labels and arrows.</li>
<li><b>Process 3D model</b>, check it in the 3D tab, and <b>Save STL</b>. Your painting and settings are saved
next to the STL, so the model can be rebuilt.</li>
</ol>
<h3>Designing for touch</h3>
<ul>
<li>Fingers find <b>edges and steps</b> far more easily than gentle slopes. Features need to be about
0.5-1 mm high and about 3 mm apart to be told apart.</li>
<li>Small round bumps can be mistaken for <b>braille dots</b>: remove or lower compact sources unless they matter.</li>
<li>Say out loud what is <b>not data</b>: markers, textures, shapes and exaggerated heights.</li>
<li>Have a <b>braille reader check every label</b> before printing a set.</li>
</ul>
<p>Every section has a <b>?</b> button that explains its settings.</p>
"""

SECTIONS = {
    "Templates": """
<p><b>Plaque format</b> sets the physical size: how long the image is, how tall the relief is, the base, the
border and the label strip. <b>Image recipe</b> sets how the image is prepared: cleaning, star removal and the
stretch. Choosing one resets just its own settings, so the two can be mixed freely.</p>
<p><b>Save template</b> stores all your current settings (but not your painting, markers or labels) to reuse
on another image; <b>Load</b> applies one. Templates also work on the command line
(<code>--plaque</code>, <code>--recipe</code>, <code>--template</code>).</p>""",

    "Model size and printer": """
<p><b>Longest side</b> is the image's size on the plaque; the border and label strip add to it.
<b>Maximum relief</b> is how tall the brightest feature stands above the base: 8-12 mm is easy to feel on a
large plaque, 5 mm on a small one. The <b>base</b> keeps the plaque rigid.</p>
<p>The <b>print bed</b> size is only used to warn you, below the settings, when the plaque will not fit
(even rotated).</p>""",

    "Image and stretch": """
<p>The <b>stretch</b> decides how brightness becomes height. Astronomical images span a huge range of
brightness, so a straight (linear) mapping usually shows only the brightest spot.</p>
<ul>
<li><b>linear</b>: height in proportion to brightness. Good for images already stretched for display.</li>
<li><b>log</b> and <b>asinh</b>: compress bright parts so faint ones get height too. asinh is gentler on the
faintest levels. For data with tiny values (e.g. Jy/beam) set <b>Soften</b>, or use noise cleaning, or they
behave almost linearly.</li>
<li><b>tactile</b>: designed for touch. See the <i>Tactile stretch</i> section.</li>
</ul>
<p><b>The chart</b> plots every part of the image from faintest (left) to brightest (right) against the height
it gets. A long flat start means most of the image is flat base; a steep end means the brightest parts take
most of the relief. Points above the line were lifted by the contrast boost.</p>
<p><b>Clip</b> ignores the faintest and brightest few percent (auto: 1%, none with cleaning or tactile).
<b>Smoothing</b> takes the grit off in model pixels. <b>Downsampling</b> averages pixels together;
the <b>smallest model pixel</b> stops the model from being finer than a printer can print, which only makes
huge files. <b>Ridge preservation</b> keeps thin filaments continuous when downsampling (cleaned maps only).</p>""",

    "Tactile stretch": """
<p>The tactile stretch turns brightness into height in a way fingers can read:</p>
<ol>
<li><b>Edge</b>: everything fainter than this is exactly flat base, so empty sky is smooth. Auto: the sky
level plus 3 times the noise (with cleaning: the island edge).</li>
<li><b>Outline step</b>: the object's edge rises straight to this fraction of the relief (default 20%), so a
finger finds the outline at once instead of a slope that fades into nothing.</li>
<li><b>Log scale</b> from the edge up to the peak, so faint outer parts and the bright centre both get room.</li>
<li><b>Peak cap</b>: the brightest 0.1% are capped, so a bright core becomes a small summit instead of taking
most of the height and leaving discs and arms flat.</li>
<li><b>Contrast boost</b>: each part is compared with its surroundings (over about 2.5% of the image) and the
difference is exaggerated. Arms, bars, knots and shells stand up from the galaxy around them; smooth shapes
barely change. 0 turns it off; 3-5 for faint arms or thin dust lanes.</li>
</ol>""",

    "Noise cleaning (radio maps)": """
<p>Radio images are full of noise. Printed as is, it feels like sandpaper, and noise peaks become bumps that
can't be told from real galaxies. Cleaning finds real emission and makes everything else perfectly flat.</p>
<ul>
<li>An <b>island</b> of emission is kept if its edge is above the <b>island edge</b> (3 sigma = 3 times the
noise) and it peaks above the <b>island peak</b> (5 sigma) somewhere, and it is at least the
<b>smallest island</b> (half a beam).</li>
<li>Emission is also searched for in a smoothed copy, where faint, spread-out lobes rise above their noise.
Add a broader <b>detection scale</b> (e.g. <code>3 9</code>) for very faint, very extended emission.</li>
<li><b>Local noise boxes</b> measure the noise piece by piece, for mosaics whose noise changes across the field.</li>
<li>The <b>beam</b> (the telescope's resolution) is read from the file; enter it if the file lacks it.</li>
</ul>
<p>The result is in units of the noise, so use the <b>asinh</b> or <b>tactile</b> stretch afterwards.</p>""",

    "Stars and compact sources": """
<p>In optical images, foreground stars of our own galaxy print as sharp spikes that feel just like
star-forming knots. Removed areas are filled smoothly from their surroundings, so no bump or dip is left.</p>
<ul>
<li><b>Gaia</b> removes known foreground stars (needs the network the first time; the list is cached). Only
stars that move, or are bright, are removed, so clusters inside a galaxy stay.</li>
<li><b>Compact sources</b> removes small sharp things Gaia misses: faint stars, globular clusters, small
background galaxies. Skip it for spirals and irregulars whose knots matter. Widen the <b>median filter</b>
for images with broad stars (about 21 for JWST pictures).</li>
</ul>
<p>Turn on <b>Show removed</b> above the 2D view to see what was taken out: cyan for Gaia stars, magenta for
compact sources. Use the <b>Protect</b> brush on anything that must stay, such as the nucleus.</p>""",

    "Paint-out, texture and markers": """
<ul>
<li><b>Paint-out</b>: <i>fill from surroundings</i> removes the painted area as if it was never there;
<i>flatten to base</i> makes it flat.</li>
<li><b>Textures</b> (dots, lines, grid) let regions be told apart by touch. Keep them clearly different from
braille: dots well over 2.5 mm apart.</li>
<li><b>Remove islands under</b> a size clears small isolated bumps that could be mistaken for braille dots
(needs a flat background: cleaning or the tactile stretch). Or keep them as low pips.</li>
<li><b>Markers</b> are smooth domes for things too faint to feel but central to the story, such as a nebula's
central star. A marker is not data: say so.</li>
</ul>""",

    "Shapes and labels on the image": """
<ul>
<li><b>Labels</b> on the image sit on a flat pad, because braille can't be read over bumps. Put them in clear
sky: a pad over a source flattens part of it (you'll be warned). Braille can't be shrunk, so keep labels short.</li>
<li><b>Lines and arrows</b> point from a label to a feature; <b>rings</b> mark an area. They follow the relief.
Raised lines read best at least 1 mm wide and 0.5-1 mm high; on coarse images, make them wider.</li>
</ul>""",

    "Border and label": """
<p>A raised <b>border</b> protects the edges and tells a finger where the image ends.</p>
<p>The <b>label strip</b> is a flat area below the image for the label: allow about 10 mm per braille line.
<b>Braille</b> is Unified English Braille, uncontracted, with standard dots (0.8 mm high). <b>Text</b> is
raised print for sighted helpers and low-vision visitors. Have a braille reader check every label.</p>""",
}

TOOLS = {
    "pan": "Pan: drag to move the view. Scroll to zoom with any tool; middle-drag also pans.",
    "paint": "Paint out: paint over anything to remove it. Set fill or flatten in 'Paint-out, texture and markers'. "
             "Ctrl/Cmd+Z undoes.",
    "protect": "Protect: paint over things star and source removal must keep, such as a galaxy nucleus.",
    "texture": "Texture: paint a raised pattern. Choose dots, lines or grid to the right; spacing and height in "
               "'Paint-out, texture and markers'.",
    "erase": "Erase: removes paint-out, protect and texture under the brush.",
    "marker": "Marker: click to place a dome; right-click one to remove it. Markers are not data: say so.",
    "label": "Label: click in clear sky to place a label; double-click a label to edit it; right-click to delete.",
    "line": "Line: click to add points; double-click or Enter to finish, Esc to cancel. Right-click a line to delete it.",
    "arrow": "Arrow: like Line, with an arrowhead at the last point. Start at the label, end at the feature.",
    "circle": "Circle: drag from the centre to draw a ring. Right-click a ring to delete it.",
}
