function [group_names, group_members] = build_groups(mkt)
% Retrieve the field names of the market structure and
% build the groups based on the sector/factor fields
%
% INPUTS:
% mkt:              Market structure
% 
% OUTPUTS:
% group_names:      Names of the groups
% group_members:    Members of the groups

% Retrieve the field names of mkt
fields_mkt = fieldnames(mkt);

% Combine the field names of mkt.sector with the second field of mkt
group_names = [fieldnames(mkt.sector); fields_mkt(2)];

sector_fields = fieldnames(mkt.sector);

for i = 1:length(group_names)
    % If the field is sector, then the members are the fields of the sector
    if i < length(group_names)
        group_members{i} = mkt.sector.(sector_fields{i});
    else
        % If the field is not sector, then the members are the values of the field
        group_members{i} = mkt.(group_names{i});
    end
end

end
