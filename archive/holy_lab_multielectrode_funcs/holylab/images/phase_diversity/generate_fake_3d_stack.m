function stack_out = generate_fake_3d_stack

% generates 3d stack
% very dumb and non robust way to do this. but works.

stack = zeros([256 256 64]);
stack = single(stack);


stack1 = create_sphere(stack, [10 108 20], 30, 150);

stack2 = create_cube(stack, [10 10 10], 20, 250);

stack3 = create_sphere(stack, [128 128 50], 10, 200);

stack4 = create_random_dots(stack, [200 200 50], 10, 2, 100);

stack5  = create_random_dots(stack, [100 100 30], 10, 5, 230);

stack6 = create_sphere(stack, [200 200 20], 20, 200);

stack7 = create_cube(stack, [100 10 50], 10, 180);

stack8 = create_sphere(stack, [200 200 10], 20, 200);

stack9 = create_sphere(stack, [10 200 20], 5, 500);

stack10 = create_sphere(stack, [200 10 20], 5, 500);

stack_out = stack1+stack2+stack3+stack4+stack5+stack6+stack7+stack8+stack9+stack10;

% figure; 
% for indx = 1:round(size(stack_out,3)/10):round(size(stack_out,3))
%     %colormap(gray(256));
%     image(stack_out(:,:,indx));
%     axis image; 
%     title(num2str(indx));
%     pause;
% end

end

function stack_out = create_sphere(stack_in, center, radius, gray_level)

stack_out = zeros(size(stack_in));
for indx1 = 1:size(stack_in, 1)
    for indx2 = 1:size(stack_in, 2)
        for indx3 = 1:size(stack_in, 3)
            if ((indx1-center(1))^2 + (indx2-center(2))^2 + (indx3-center(3))^2) <= radius^2
                stack_out(indx1, indx2, indx3) = stack_in(indx1, indx2, indx3) + gray_level;
            end
        end
    end
end

end

function stack_out = create_cube(stack_in, start_point, edge_size, gray_level)

stack_out = zeros(size(stack_in));
for indx1 = start_point(1):start_point(1) + edge_size
    for indx2 = start_point(2):start_point(2) + edge_size
        for indx3 = start_point(3):start_point(3) + edge_size
            stack_out(indx1, indx2, indx3) = stack_in(indx1, indx2, indx3) + gray_level;       
        end
    end
end

end


function stack_out = create_random_dots(stack_in, start_point, edge_size, spacing, gray_level)

stack_out = zeros(size(stack_in));
for indx1 = start_point(1):spacing:start_point(1) + edge_size
    for indx2 = start_point(2):spacing:start_point(2) + edge_size
        for indx3 = start_point(3):spacing:start_point(3) + edge_size
            stack_out(indx1, indx2, indx3) = stack_in(indx1, indx2, indx3) + round(rand)*gray_level;       
        end
    end
end

end


%function stack_out = create_stick_cube

%function stack_out = creat_pyramid(stack_in, 


