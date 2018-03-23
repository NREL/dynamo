function [results] = zmq_client(data, ipaddress)

        message = savejson('', data);

        persistent ctx;
        % set original servo positions to zero the
        % first time this function is invoked
        if isempty(ctx)
            disp('Creating CTX')
            import org.jeromq.ZMQ;
            import org.jeromq.ZContext;
            ctx = ZContext();
        end

        socket = ctx.createSocket(ZMQ.REQ);

        socket.connect(ipaddress);

        socket.send(unicode2native(message));

        message = socket.recv();

        m = native2unicode(message)';

        results = loadjson(m);

        socket.close();

end
