function steps = fn_unmap_nds_in_steps(steps)
%undoes the effect of fn_remap_nds_in_steps
for s = 1:numel(steps)
    if isfield(steps{s}.load, 'orig_frc_nds')
        steps{s}.load.frc_nds = steps{s}.load.orig_frc_nds;
        steps{s}.load = rmfield(steps{s}.load, 'orig_frc_nds');
    end
    if isfield(steps{s}.load, 'orig_dsp_nds')
        steps{s}.load.dsp_nds = steps{s}.load.orig_dsp_nds;
        steps{s}.load = rmfield(steps{s}.load, 'orig_dsp_nds');
    end
    if isfield(steps{s}.mon, 'orig_nds')
        steps{s}.mon.nds = steps{s}.mon.orig_nds;
        steps{s}.mon = rmfield(steps{s}.mon, 'orig_nds');
    end
end