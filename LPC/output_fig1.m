
       mean_cis=[mean_cis,mean_ci];
    std_cis=[std_cis,std_ci];
      
   for i=1:7  
       subplot(1,2,1);
       grid on;
       hold on; 
       switch(i) 
      case 1
         errorbar([1:10],  mean_cis(:,i),std_cis(:,i), 'color',[1 0 1]);
      case 2
         errorbar([1:10],  mean_cis(:,i),std_cis(:,i),  'color',[0 1 1]);
      case 3
         errorbar([1:10], mean_cis(:,i),std_cis(:,i),  'color',[1 0 0]);
      case 4
         errorbar([1:10],  mean_cis(:,i),std_cis(:,i),  'color',[0 1 0]);
      case 5
         errorbar([1:10],  mean_cis(:,i),std_cis(:,i),  'color',[0 0 1]);
      case 6
         errorbar([1:10],  mean_cis(:,i),std_cis(:,i), 'color',[0 0.5 0.3]);
      case 7
         errorbar([1:10],  mean_cis(:,i),std_cis(:,i),  'color',[0.6 0.2 0]);
       end
  subplot(1,2,2);
  grid on;
     hold on;
  switch(i)
      case 1
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[1 0 1]);
      case 2
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0 1 1]);
      case 3
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[1 0 0]);
      case 4
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0 1 0]);
      case 5
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0 0 1]);
      case 6
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0 0.5 0.3]);
      case 7
         errorbar([1:10], mean_num_shd(i,:), std_num_shd(i,:), 'color',[0.6 0.2 0]);
  end
   end