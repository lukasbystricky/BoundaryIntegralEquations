function [dirVs, neuVs] = create_bc(xq, yq, f_handle, g_handle)


dirVs = zeros(size(xq));
neuVs = zeros(size(xq));

for i = 1:length(xq)
    dirVs(i) = f_handle(xq(i),yq(i));
    neuVs(i) = g_handle(xq(i),yq(i));
end

