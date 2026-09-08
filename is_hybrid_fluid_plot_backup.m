function tf = is_hybrid_fluid_plot(out)
    tf = false;
    if ~isfield(out, 'stopStep') || out.stopStep < 1 || ...
            ~isfield(out, 'fluidHist') || numel(out.fluidHist) < out.stopStep || ...
            isempty(out.fluidHist{out.stopStep})
        return;
    end
    fluidPlot = out.fluidHist{out.stopStep};
    tf = isfield(fluidPlot, 'meshType') && strcmpi(fluidPlot.meshType, 'hybrid_gap1d_exterior2d');
end