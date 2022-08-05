function features = ssr_keep_features(features_in,keep_index)
  features = features_in;
  features.Xmin = features.Xmin(keep_index,:);
  features.Xmax = features.Xmax(keep_index,:);
  features.Xcenter = features.Xcenter(keep_index,:);
  features.offset = features.offset(:,:,keep_index);
  features.T = features.T(:,keep_index);
  features.S = features.S(keep_index);
  