function color_value = generatecolormapthreshold(threshold,colour)
    arguments
        threshold (:,:) {mustBeNumeric}
        colour (:,:) {mustBeNumeric, mustBeGreaterThanOrEqual(colour,0), mustBeLessThanOrEqual(colour,1)}
    end
    validateattributes(threshold,{'numeric'},{'size', [1 NaN]},mfilename,"Threshold",2)
    validateattributes(colour,{'numeric'},{'size', [NaN 3]},mfilename,"Colour",2)
    if length(threshold)-1 ~= size(colour,1)
        error("threshold and colour are not match. Either more colour or more threshold is assigned")   
    else
        color_value=[];
        for i = 1:1:length(threshold)-1
            createrange = threshold(i):threshold(i+1);
            for w = 1:1:length(createrange)-1
               color_value = [color_value;colour(i,:)];
            end
        end
    end
end