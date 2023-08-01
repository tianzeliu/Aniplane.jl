a = zeros(10, 1)
Threads.@threads for i = 1:10
           b = Threads.threadid()
           a[i] = b^2
end

print(a)