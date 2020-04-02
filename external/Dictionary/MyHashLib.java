public class MyHashLib {
    
    // FNV1A 64 bit
    public static long fnvHash64(byte[] s) {
        long h = 0xcbf29ce484222325L;
        for (byte b : s) {
            h ^= b;
            h *= 0x100000001B3L;
        }
        
        return h;
    }

    public static long fnvHash64(String s) {
        return fnvHash64(s.getBytes());
    }
    
    // FNV1A 32 bit
    public static int fnvHash32(byte[] s) {
        int h = 0x811C9DC5;
        for (byte b : s) {
            h ^= b;
            h *= 0x1000193L;
        }
        
        return h;
    }

    public static int fnvHash32(String s) {
        return fnvHash32(s.getBytes());
    }
    
    // Jenkins Hash
    public static int jenkinsHash(byte[] s) {

        int h = 0;
        for (byte b : s) {
            h += b;
            h += (h << 10);
            h ^= (h >>> 6);
        }
        h += (h << 3);
        h ^= (h >>> 11);
        h += (h << 15);
        
        return h;
    }
    
    public static int jenkinsHash(String s) {
        return jenkinsHash(s.getBytes());
    }
}