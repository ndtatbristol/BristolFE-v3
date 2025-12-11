function steps = fn_remap_nds_in_steps(steps, new_nds)
%Used if it is necessary to remap all node references in steps (e.g. after
%node removal or node sorting. Original node numbers are copied to orig_*
%variables for later restoration if required
for s = 1:numel(steps)
    if isfield(steps{s}.load, 'frc_nds')
        steps{s}.load.orig_frc_nds = steps{s}.load.frc_nds;
        steps{s}.load.frc_nds = fn_remap_matrix(steps{s}.load.frc_nds, new_nds);
    end
    if isfield(steps{s}.load, 'dsp_nds')
        steps{s}.load.orig_dsp_nds = steps{s}.load.dsp_nds;
        steps{s}.load.dsp_nds = fn_remap_matrix(steps{s}.load.dsp_nds, new_nds);
    end
    if isfield(steps{s}.mon, 'nds')
        steps{s}.mon.orig_nds = steps{s}.mon.nds;
        steps{s}.mon.nds = fn_remap_matrix(steps{s}.mon.nds, new_nds);
    end
end
end