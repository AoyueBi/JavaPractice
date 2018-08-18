/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package analysis.collaboration;

/**
 *
 * @author Fei Lu <fl262@cornell.edu>
 */
public class CollaborationGo {
    
    public CollaborationGo () {
        //this.vtePipe();
        //this.jianbingPipe();
        //this.joyEiPipe();
        //this.w22Pipe();
        //this.b104Pipe();
        //new Randy();
        //new Feng();
        //new Tao();
        new Philip();
    }
    
    public void b104Pipe () {
        new B104();
    }
    
    public void w22Pipe () {
        new W22();
    }
    
    public void joyEiPipe () {
        new JoyEi();
    }
    
    public void jianbingPipe () {
        new Jianbing();
    }
    
    public void vtePipe () {
        new VTE();
    }
    
    public static void main (String[] args) {
        new CollaborationGo();
    }
}
