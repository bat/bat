inkscape = inkscape --export-area-drawing

%.png: %.svg
	$(inkscape) --export-png=$@ --export-width 500 $^

histogram.png: histogram.svg
	$(inkscape) --export-png=$@ --export-width 250 $^

%.pdf: %.svg
	$(inkscape) --export-pdf=$@ $^
