function save_pdf_figure(fig, file_name)
% Save figures as tightly bounded vector PDFs for inclusion in LaTeX.

if nargin < 1 || isempty(fig)
    fig = gcf;
end

set(fig, 'Color', 'w');

try
    exportgraphics(fig, file_name, 'ContentType', 'vector', 'BackgroundColor', 'white');
catch
    set(fig, 'Units', 'inches');
    pos = get(fig, 'Position');
    set(fig, 'PaperUnits', 'inches');
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'PaperPosition', [0 0 pos(3) pos(4)]);
    set(fig, 'PaperSize', [pos(3) pos(4)]);
    print(fig, '-dpdf', '-painters', file_name);
end
end
